# -*- coding: cp1252 -*-
# transduction.py
import numeric_id
import collections, re, pdb, copy, datetime, cPickle
from collections import defaultdict

"""
    Takes an English Sig string (directions to patients on how to administer a medicine) and tries to break down it's meaning into components via a process  called "transduction".

    The main procedure is  transduce_sig(english_sig_string).

    The algorithm consists of 4 main parts:

    apply_struc_identification_rules(sig):  Identifies linear subsegments of the string as meaningful chunks (Strucs). There are 3 dozen such Strucs, e.g. Dose, Freq, Directive, etc.
                                            The Strucs and their relationship to Sems are documented in Transduction/struc_types.xlsx.

    create_sems(sig):   Takes the sequence of Strucs and tries to create from them a semantic tree (Sem) that represents the scheduling during the day and between days of basic directions.
                        The basic (atomic) administation structure is called an AdminEvent which represents direction to take/use some drug at one time. Several AdminEvents that represent
                        multiple administrations in one day are grouped into a Schedule. Multiple schedules may represent multiple drug administration regimes, and are grouped into Instructions.
                        Instruction/Schedule/Event
                        For example, "take 2 tabs now, then for the next 2 weeks take 1 tab every morning and 2 tabs before bedtime" is represented as 3 events: 1 event in 1 schedule and 2 events in
                        the second schedule.

    modify_sems(sig):   Reasons about the Sems, trying to figure out and modify such things as the main directive, the main drug form, etc. This is needed because for many if not most sigs a lot of information
                        is missing or is implicit.

    break_sig_into_atoms(sig):  Final step where the Sems from the Sig are mapped into semantic represntation (Sems) of the dictionary entries (called atoms). The resulting matches are grouped into Versions and
    stored in sig.versions. The top several versions can be instantiated for specific locales from the pickled dictionary at code/pickles/atoms_dictionary.pkl.


    To debug:
    There are several quasi-unit tests available:
    test_numeric_id.py              Tests numeric identification. numeric_id.py is the only local imported module for transduction.
    test_struc_identification.py:   tests apply_struc_identification_rules().
    test_sems.py:                   Tests create_sems()
    test_transduction.py            The most useful and important module. It runs the algorithm against the >10K sigs and their existing (or prior) transductions localed in code/gold_standard_pairs.xlsx.
                                    It is not so much that the transductions recorded in that excel file are perfect or desirable. Rather, because they include all popular sigs and examples of every
                                    feature we have seen and dealth with so far, and because of their sheer number, we expect that any non-trivial modification to the code will result in some changed transduction.
                                    So this procedure allows us to see what is the effect of any change we made on past transductions.
                                    As any new features are added to the code, make sure that at least a dozen or so different raw sigs with that feature are added to the spreadsheet.

    Tips for debugging:
    -   A useful way to understand how a particular sig got transduced or matched to atoms in the way it did is to use
        sigs = test(raw_sigs = raw_sigs[0:1], debugging_flag = True, details_to_show = 0, limit_to_incompletely_parsed = False, return_sigs = True, print_developer_msgs = False)
        Then focus on sig = sigs[0] and run
            for struc in sig.parses[0].strucs: print struc
        Then for each struc printed look at the rules_used property list.

    -   To see how atoms are represented, use show_atom(raw_atom_string).


    Maintenance:
    -   When making any changes to the code, always run test_transduction.py to see it's effect on the prior transductions. When adding any features, make sure that at least a few sigs exhibiting these features
        are added to the code/gold_standard_pairs.xlsx spreasheet (and thereby are in the future used by test_transduction.py).
    -

    -   Relevant files:
        Used by code:
            code/pickles/atoms_dictionary.pkl
            code/pickles/tentative_atomic_dictionary.pkl
            code/numeric_id.py

        Used to process new dictionary entries and atoms:
            - Translations/dictionary_tentative_atomic.xls
            - code/dictionary.py

        Tests:
            test_numeric_id.py
            test_struc_identification.py
            test_sems.py
            test_transduction.py

        Data:
            Transduction/struc_types.xlsx
            code/gold_standard_pairs.xlsx
            sig_codes/SIG_Abbreviations.xlsx

        Debugging and building datasets:
            code/experiment_with_transduction.py
                Use export_transductions_quick() and get_patterns() to produce transductions of either a large set of sigs or to zoom in on specific types of sigs.
            code/build_datasets.py
                Use to clean up raw query logs, merge corpora to create frequency sorted raw sigs, and to run transduction on corpora (use transduce_sigs_with_freqs())


"""


##### CONSTANTS #######

debug = False       # transduce() can set debugging on with parameter debugging_flag = True
RULES = []          # collection of all Rules that may be applied to a sig

[arm, ar, bn, ch, de, en, fr, gr, hi, ht, it] = ["hy_AM", "ar_AR", "bn_BD", "zh_CN", "de_DE", "en_US", "fr_FR", "el_GR", "hi_IN", "ht_HT", "it_IT"]
[jp, kr, pl, pt, ru, sp, tg, tr, tw, vt] = ["ja_JP", "ko_KR", "pl_PL", "pt_BR", "ru_RU", "es_US", "tl_PH", "tr_TR", "zh_TW", "vi_VN"]
[amh, bur, far, nep, pas, rom, som, swa] = ["am_ET", "my_MM", "fa_IR", "ne_NP", "ps_AF", "ro_RO", "so_SO", "sw_KE"]
available_locales = [arm, ar, de, en, fr, gr, it, ru, sp, ch, ht, kr, pl, pt, tg, tw, bn, hi, vt, amh, bur, far, nep, pas, rom, som, swa]
locales_using_period_for_decimal = set([en, sp, tg, kr, nep, bur])

canonical_struc_order = ('THEN_CHRONO', 'AND_CONJ', 'TAPER', 'DIRECTIVE', 'ANAPHORA', 'DOSE', 'SUBSTRATE', 'ROUTE', 'SITE', 'VEHICLE', 'DIR_ATTRIBUTE',
                         'TIMING', 'SPECIFIC_DAY', 'FREQ', 'PERIODICITY', 'CALENDAR_EVENT', 'DURATION',
                         'DISCARD_REMAINDER', 'STOP_CONDITION', 'AS_NEEDED', 'REPEAT', 'MAXDOSE', 'AS_DIRECTED', 'SEE_INSTRUCTIONS', 'INDICATION', 'MISCELLANEOUS')

errors_list = []                    # List of DeveloperMessages objects describing errors/msgs for the current sig being persed
collective_errors_all_sigs = []     # List of all errors_lists for the various sigs

##### CLASSES ########

class UnparseableSig(Exception):
    def __init__(self, value, raw_sig):
        self.message = value
        self.raw_sig = raw_sig
    def __str__(self):
        txt = ('Error: |%s|  In raw sig: |%s|' % (self.message, self.raw_sig))
        return txt

class DeveloperMessages(object):
    def __init__(self, msg_string, msg_type, error_level = None):
        self.msg = msg_string
        self.msg_type = msg_type
        self.error_level = error_level if error_level else 0

    @classmethod
    def print_developer_messages(cls):
        if collective_errors_all_sigs:
            print('\n------\nDEVELOPER MESSAGES: Found %d sigs with errors' % len(collective_errors_all_sigs))
            for sig_num, errors_list in enumerate(collective_errors_all_sigs):
                raw_sig = errors_list[0].msg
                print('Sig %d. For raw sig %s found %d messages' % (sig_num, raw_sig, len(errors_list) - 1))
                for i, error_msg in enumerate(errors_list[1:]):
                    print('Message %d of type %s, error level %d:\n      %s' % (i, error_msg.msg_type, error_msg.error_level, error_msg.msg))


class Sig(object):
    """
    The sig object.

    raw_sig         String: actual raw txt as submitted by the customer
    cleaned_raw_sig String: after light string cleansing and normalization.
    dictionary      Dict(): is used for atoms from the dictionary. Is None for all customer sigs. For atoms, it is a dictionary mapping locale to the string translating this atom for this locale.
    parses          List of Parse objects, which are the currently considered partial passes.
    versions        List of Version objects, which contain the various maximal ways the raw_sig can be matched to lists of Atoms from the dictionary.
    """

    def __init__(self, raw_sig, cleaned_raw_sig):

        self.raw_sig = raw_sig
        self.dictionary = None
        self.cleaned_raw_sig = cleaned_raw_sig
        self.parses = [Parse(self)]
        self.versions = []

    def show(self, details_to_show = None):
        """
        details_to_show     integer parameter -- higher number implies more details. Currently 3 max
        """

        if not details_to_show:
            details_to_show = 4

        rows = []
        rows.append('Raw Sig:       %s' % self.raw_sig)
        if details_to_show > 1:
            rows.append('Cleaned Sig:   %s' % self.cleaned_raw_sig)
        rows.append('Number of Parses: %d' % len(self.parses))
        for parse_num, parse in enumerate(self.parses):
            parse_descr = ('---> Parse %d with %d Instructions and %d Strucs. Cycles run: %d (%s; %s). Max recursive calls for atom match: %d' %
                           (parse_num, len(parse.instructions), len(parse.strucs), parse.number_of_cycles_of_rule_application,
                            'fully parsed to strucs' if parse.is_fully_parsed() else 'incomplete struc parse',
                            'fully segmented to sems' if parse.is_fully_sem_segmented() else 'incomplete sem segmentation',
                             parse.debug_max_recursive_calls))
            labels = ('\n     Labels: %s' % parse.flatten())
            sem_descr_basic = ('\n     SEMS:   ' + parse.show_struc_assignment_to_sems(include_struc_details = False, omit_coords = False) +
                               '\n     SEMS:   ' + parse.show_struc_assignment_to_sems(include_struc_details = True, omit_coords = False) +
                               '\n     Flattened SEM:   ' + parse.flatten_sems() +
                               '\n     Primary Directive: >' + (parse.instructions[0].primary_directive.value if parse.instructions[0].primary_directive else '') +
                                       '< Primary Form: >' + (parse.instructions[0].primary_form.value if parse.instructions[0].primary_form else '') +
                                       '< Primary Route: >' + (parse.instructions[0].primary_route.value if parse.instructions[0].primary_route else '') + '<' +
                               ('\n     Key set (%d atoms matched):    %s' %
                                (len(get_sorted_compatible_atoms(parse)), get_normalized_key_set(parse.strucs)))
                               )

            sem_descr = labels + sem_descr_basic
            if details_to_show > 1:
                sem_descr += '\n\n' + parse.sem_scope_data.pprint() + '\n'
                sem_descr += padlines(parse.show_sems(), 5)
            rules_descr = ('\n     Rules: %d: %s' % (len(parse.rules_utilized), '; '.join(parse.rules_utilized)))
            strucs = ('\n     %d Strucs:%s' % (len(parse.strucs), padlines(parse.__str__(), 10)))
            row = parse_descr + sem_descr
            if details_to_show > 2:
                row += rules_descr

            if details_to_show > 3:
                row += strucs
            #if details_to_show > 1:
            #    row += '\n' + labels + ('\n     Labels: %s' % parse.pprint_strucs()) + sem_descr_basic
            if details_to_show > 1:
                row += '\n\n----------- START VERSIONS: Found %d versions %s. ' % (len(self.versions), ' but will show 2' if len(self.versions) > 2 else '')
                row += ('Fully parsed to strucs' if parse.is_fully_parsed() else 'incomplete struc parse') + "|"
                row += ('fully segmented to sems' if parse.is_fully_sem_segmented() else 'incomplete sem segmentation')
                row += padlines('\nRaw Sig: -->%s<--' % self.raw_sig, 20)
                for ver_num, version in enumerate(self.versions[:4]):
                    row += padlines('\n' + version.show(ver_num), 5)
                row += '\n----------- END VERSIONS\n'
            rows.append(row)

            #rows.append('---> Parse %d (%s): \n     Rules: %d: %s\n     Labels: %s\n     %d Strucs:%s' %
            #            (parse_num, 'fully parsed' if parse.is_fully_parsed() else 'incomplete parse', len(parse.rules_utilized), '; '.join(parse.rules_utilized),
            #             parse.flatten(), len(parse.strucs), padlines(parse.__str__(), 10)))
        return '\n' + '\n'.join(rows)

    def pprint(self, details_to_show = None):
        print self.show(details_to_show = details_to_show)

    def __str__(self):
        return self.show(details_to_show = None)

class Parse(object):
    """
    A partial state of processing of the sig.

    strucs          LIST: A list of constituent structures at the current parsing stage. The strucs list may
                    look like this: DIR NUM tablet  (where each capitalized word represents a STRUC object
                    with the capitalized word as its label).

    instructions    LIST: A list of Instructions, each representing a separate thought, e.g. administration schedule (most typical),
                    request to make an appointment, etc.

    changed_on_last_pass                    True/False: A flag used to signal whether the last pass through a set of rules modified the parse or not.
    number_of_cycles_of_rule_application    INT. For debugging
    debug_max_recursive_calls               INT. For debugging: counts the maximal number of calls to find_maximal_instatom_groups_recursively() for an event.
    """

    def __init__(self, sig, strucs = None):
        """
        """
        self.sig = sig
        self.strucs = strucs if strucs else self.initialize_strucs()
        self.instructions = []
        self.rules_utilized = []        #   List of Rules that made any changes to strucs or any struc in strucs.
        self.changed_on_last_pass = True
        self.number_of_cycles_of_rule_application = 0
        self.debug_max_recursive_calls = 0

    def initialize_strucs(self):
        txt = self.sig.cleaned_raw_sig
        strucs = [Struc(txt)]
        return strucs

    def flatten(self, start_index = None):
        if not start_index:
            start_index = 0
        result = get_struc_labels(self.strucs[start_index:], delimiter = '', omit_spaces_and_punctuation = False)

        return result

    def copy(self):
        """ Returns a copy of self, made as shallowly as possible."""

        new_strucs = [copy.copy(struc) for struc in self.strucs]
        new_parse = Parse(self.sig, new_strucs)
        new_instructions = [copy.copy(instruction) for instruction in self.instructions]
        new_parse.instructions = new_instructions
        new_parse.rules_utilized = copy.copy(self.rules_utilized)
        return new_parse

    def is_fully_parsed(self):

        for struc in self.strucs:
            if struc.is_unparsed_non_punctuation_txt():
                return False
        return True

        #pattern_non_text = re.compile('[^a-z0-9]')
        #unparsed_txt = ''.join([struc.label for struc in self.strucs if struc.is_unstructured_text() and not struc.label.isspace()])
        #remainder = pattern_non_text.sub('', unparsed_txt)
        #remainder = unparsed_txt.replace('.', '').replace(',', '')
        #if not remainder or remainder.isspace():
        #    return True
        #else:
        #    return False

    def position2struc(self, position):
        """ Returns the struc whose label's position in the self.flatten() string covers the position param.

            E.g., if self.flatten() = 'DIR QUANT FORM ...', and position is between 4 and 8, then the returned
            value is the struc corresponding to QUANT.
        """

        end_pos = -1
        for struc_num, struc in enumerate(self.strucs):
            start_pos = end_pos + 1
            end_pos += len(struc.label)
            if position >= start_pos and position <= end_pos:
                return struc
        return None

    def pointwise_equal_by_strucs(self, another_parse):
        """ Two parses are pointwise equal if their self.strucs lists are pointwise equal, struc by struc """

        if len(self.strucs) != len(another_parse.strucs):
            return False
        for i, struc in enumerate(self.strucs):
            a_struc = another_parse.strucs[i]
            if not struc.pointwise_equal(a_struc):
                return False
        return True

    def resegment(self, new_strucs, start_pos_new_strucs, end_pos_new_strucs, rule_name = None):
        """ Resegments parse.strucs list by removing everything in parse.flatten() from start_pos_new_strucs to end_pos_new_strucs
            and replacing these text fragments and strucs with a list of new_strucs.

            Note that end_pos_new_strucs is the actual position of the last character in parse.flatten(), i.e.
            new_strucs replace parse.flatten()[start_pos_new_strucs: end_pos_new_strucs + 1].

            If start_pos_new_strucs or end_pos_new_strucs falls in the middle of a struc, that struc has to be a pure text structure.
        """

        # First, determine how the new_strucs fit with the existing strucs, i.e. where do they split the existing strucs and how.
        end_pos = -1
        for struc_num, struc in enumerate(self.strucs):
            start_pos = end_pos + 1
            end_pos += len(struc.label)

            # We need to find after which struc in the self.strucs list to start inserting new_strucs, and
            # at which struc in the self.strucs the inserting of new_strucs ends. I.e. where is the gap in the
            # strucs list into which to insert new segments.
            # These start and end strucs we call start_gap_struc and  end_gap_struc.
            # We also need to know if start_gap_struc and end_gap_strucs get split up by the new_strucs, and
            # if so, at what characters.

            if start_pos_new_strucs >= start_pos and start_pos_new_strucs <= end_pos:
                start_gap_struc_num = struc_num
                start_gap_struc = struc
                start_gap_start_char_num = start_pos_new_strucs - start_pos         # start_gap_start_char_num is the position of the split
                                                                                    # within the start_gap_struc (which has to be text if
                                                                                    # start_gap_start_char_num > 0)
                start_gap_struc_split_up_flag =  start_gap_start_char_num > 0
                if start_gap_struc_split_up_flag and not struc.is_unstructured_text():
                    raise Exception('Bad resegmentation request. Trying to split up segment #%d with label %s at char %d' %
                                    (struc_num, struc.label, start_gap_start_char_num))

            if end_pos_new_strucs >= start_pos and end_pos_new_strucs <= end_pos:
                end_gap_struc_num = struc_num
                end_gap_struc = struc
                end_gap_start_char_num = end_pos_new_strucs - start_pos + 1         # end_gap_start_char_num is the position of the first char after the split
                                                                                    # within the end_gap_struc (which has to be text if
                                                                                    # end_gap_start_char_num < end_pos)
                end_gap_struc_split_up_flag =  end_pos_new_strucs < end_pos
                if end_gap_struc_split_up_flag and not struc.is_unstructured_text():
                    raise Exception('Bad resegmentation request. Trying to split up segment #%d with label %s at char %d' %
                                    (struc_num, struc.label, end_gap_start_char_num))


        # Second, modify existing strucs as needed, create new ones out of existing ones if needed, delete replaced strucs, and
        # insert the new strucs into the sequence parse.strucs.
        if start_gap_struc is end_gap_struc:
            # We are splitting up a single text struc and inserting the new_strucs in the middle of it.
            old_label = start_gap_struc.label
            if start_gap_struc_split_up_flag and end_gap_struc_split_up_flag:
                # We need to create a new tail struc which will have the tail end of the text of the struc we splitting up.
                start_gap_struc.label = old_label[:start_gap_start_char_num]
                new_tail_struc = Struc(label = old_label[end_gap_start_char_num:])
                new_strucs.append(new_tail_struc)
                self.strucs[start_gap_struc_num + 1: start_gap_struc_num + 1] = new_strucs
            elif start_gap_struc_split_up_flag:
                # Insert the new_strucs after the struc we are splitting
                start_gap_struc.label = old_label[:start_gap_start_char_num]
                self.strucs[start_gap_struc_num + 1: start_gap_struc_num + 1] = new_strucs
            elif end_gap_struc_split_up_flag:
                # Insert the new_strucs before the struc we are splitting
                end_gap_struc.label = old_label[end_gap_start_char_num:]
                self.strucs[start_gap_struc_num: start_gap_struc_num] = new_strucs
            else:
                # We are replacing the struc we are splitting
                self.strucs[start_gap_struc_num: start_gap_struc_num + 1] = new_strucs
        else:
            # If you are splitting existing strucs in the middle, modify their labels.
            if start_gap_struc_split_up_flag:
                start_gap_struc.label = start_gap_struc.label[:start_gap_start_char_num]
            if end_gap_struc_split_up_flag:
                end_gap_struc.label = end_gap_struc.label[end_gap_start_char_num:]

            if start_gap_struc_split_up_flag and end_gap_struc_split_up_flag:
                # Keep both start_gap and end_gap strucs and insert the new_strucs in between
                self.strucs[start_gap_struc_num + 1: end_gap_struc_num] = new_strucs
            elif start_gap_struc_split_up_flag:
                # We are subsuming the end_gap_struc, so delete it.
                self.strucs[start_gap_struc_num + 1: end_gap_struc_num + 1] = new_strucs
            elif end_gap_struc_split_up_flag:
                # We are subsuming the start_gap_struc, so delete it.
                self.strucs[start_gap_struc_num: end_gap_struc_num] = new_strucs
            else:
                # Delete both start_gap_struc and end_gap_struc
                self.strucs[start_gap_struc_num: end_gap_struc_num + 1] = new_strucs
        if rule_name:
            for new_struc in new_strucs:
                new_struc.rules_used.append(rule_name)

    def add_new_DrugAdmin(self):
        """ Adds a new empty DrugAdmin to the parse.instructions and returns the new instruction  """

        instruction = DrugAdmin(self)
        self.instructions.append(instruction)

        return instruction

    def get_sem_componenets(self):
        """ Returns a list of all strucs that are properties of all Instructions, Schedules and Events that are in the parse.
            First come direct properties of the events in the first schedule of the first instruction, then properties of the first schedule,
            then properties of events of the second schedule, then properties of the second schedule,
            etc, and finally direct properties of the first instruction. Then come the same for the second instruction, etc.
        """

        components = []
        for instruction in self.instructions:
            components.extend(instruction.get_recursive_componenets())

        return components


    def validate_full_sem_segmentation(self):
        """ Validate that sem segments indeed partition all of parse.strucs into non-intersecting instructions, schedules, admin_events
            and that these are numbered properly in non-decreasing order.

            Returns the list of errors as strings

            Used primarily in test_sems.py module.
        """

        strucs = self.strucs
        errors = []
        last_admin_event_pos = 0
        last_sched_pos = 0
        last_inst_pos = 0
        sched_pos = None
        admin_event_pos = None

        for i, struc in enumerate(strucs):
            sem = struc.accounted_for_by_sem
            if struc.is_unparsed_non_punctuation_txt():
                #error = 'Struc #%d (%s: %s) is unlabeled text' % (i, struc.label, struc)
                #errors.append(error)
                continue
            elif struc.is_space_or_punctuation_only():
                continue
            if not sem:
                rules_used = ''.join(struc.rules_used)
                if 'removed_from_sem' not in rules_used and  'assigned_to_dose' not in rules_used:
                    error = 'Struc #%d (%s: %s) is not accounted by a sem' % (i, struc.label, struc)
                    errors.append(error)
                continue
            coords = sem.get_coordinate_in_sem()    # returns a tuple of length 1 to 3. E.g. for AE: (n, m, k) where
                                                    # n is pos of AE in sched.events, M position of sched in instructions, k pos of inst in parse.instructions
            inst_pos = coords[-1]
            if inst_pos < last_inst_pos:
                error = ('Struc #%d (%s: %s) belongs to instruc #%d, though previously saw inst #%d' %
                         (i, struc.label, struc, inst_pos, last_inst_pos))
                errors.append(error)
            elif inst_pos == last_inst_pos + 1:
                # we just started a new instruction, so reset the last_schedule and last_admin_event counters
                last_sched_pos = sched_pos = 0
                last_admin_event_pos = admin_event_pos = 0
            elif inst_pos > last_inst_pos + 1:
                error = ('Struc #%d (%s: %s) belongs to instruc #%d, though previously saw inst more than 1 removed: #%d' %
                         (i, struc.label, struc, inst_pos, last_inst_pos))
                errors.append(error)

            if len(coords) > 1:
                sched_pos = coords[-2]
                if inst_pos > last_inst_pos and sched_pos != 0:
                    # We started a new instruction and this is the first schedule we saw of that inst. Better be schedule #0
                    error = ('Struc #%d (%s: %s) belongs to instruc #%d/schedule #%d, though we just started a new instruction' %
                             (i, struc.label, struc, inst_pos, sched_pos))
                    errors.append(error)
                elif inst_pos == last_inst_pos and sched_pos < last_sched_pos:
                    # We are on the same instruction, so schedule positions should be non-decreasing
                    error = ('Struc #%d (%s: %s) belongs to instruc #%d/schedule #%d, though previously saw schedule #%d' %
                             (i, struc.label, struc, inst_pos, sched_pos, last_sched_pos))
                    errors.append(error)
                elif inst_pos == last_inst_pos and sched_pos == last_sched_pos + 1:
                    # we just started a new schedule, so reset the last_admin_event counter
                    last_admin_event_pos = admin_event_pos = 0
                elif sched_pos > last_sched_pos + 1:
                    error = ('Struc #%d (%s: %s) belongs to instruc #%d/schedule #%d, though previously saw schedule more than 1 removed: #%d' %
                             (i, struc.label, struc, inst_pos, sched_pos, last_sched_pos))
                    errors.append(error)

            if len(coords) > 2:
                admin_event_pos = coords[-3]
                if sched_pos > last_sched_pos and admin_event_pos != 0:
                    # We started a new schedule and this is the first admin_event we saw of that inst. Better be schedule #0
                    error = ('Struc #%d (%s: %s) belongs to instruc #%d/schedule #%d/admin_event #%d though we just started a new schedule' %
                             (i, struc.label, struc, inst_pos, sched_pos, admin_event_pos))
                    errors.append(error)
                elif sched_pos == last_sched_pos and admin_event_pos < last_admin_event_pos:
                    # We are on the same schedule, so admin_event positions should be non-decreasing
                    error = ('Struc #%d (%s: %s) belongs to instruc #%d/schedule #%d/admin_event #%d, though previously saw admin_event #%d' %
                             (i, struc.label, struc, inst_pos, sched_pos, admin_event_pos, last_admin_event_pos))
                    errors.append(error)
                elif admin_event_pos > last_admin_event_pos + 1:
                    error = ('Struc #%d (%s: %s) belongs to instruc #%d/schedule #%d/admin_event #%d, though previously saw admin_event more than 1 removed: #%d' %
                             (i, struc.label, struc, inst_pos, sched_pos, admin_event_pos, last_admin_event_pos))
                    errors.append(error)

            last_inst_pos = inst_pos
            if sched_pos is not None:
                last_sched_pos = sched_pos
            if admin_event_pos is not None:
                last_admin_event_pos = admin_event_pos

        return errors


    def is_fully_sem_segmented(self):
        errors = self.validate_full_sem_segmentation()
        if errors:
            return False
        else:
            return True

    def show_struc_assignment_to_sems(self, include_struc_details = None, omit_coords = None):
        """ Returns a string of struc labels with added suffix of the type '|XYZ_N_M_K|' where
            XYZ is AE or SCH or INSTDA if the struc is assigned to, correspondingly, AdminEvent, Schedule, or Instruction (of type DrugAdmin) itself,
            and (N, M, K) are the coordinates of the corresponding sem structure in Events/Schedules/Instructions hierarchy.

            If struc is unassigned to any part of an Instruction, the struc is represented as label + |N|

            E.g., typical result maybe 'DIRECTIVE|AE_0_1_0| unparsedblah DOSE|N|'

            omit_coords         T/F. If True, instead of |XYZ_N_M_K| returns |Y|

            include_struc_details   T/F. If True, includes struc properties along with each struc label (e.g. instead of 'DIRECTIVE|AE_0_1_0|'
                                    it would return 'DIRECTIVE|take-as_single_dose|AE_0_1_0|
        """
        if include_struc_details:
            txt = ''.join([struc.label if struc.is_unstructured_text() else  struc.label + '|' + struc.quick_print_struc() + '{' +
                           struc.show_accounted_for_by_sem(omit_coords = omit_coords) + '}' for struc in self.strucs])
        else:
            txt = ''.join([struc.label + '{' + struc.show_accounted_for_by_sem(omit_coords = omit_coords) + '}' for struc in self.strucs])
        return txt

    def show_sems(self):
        return '\n' + '\n'.join(['--INSTRUCTION #%d: %s' % (inst_num, padlines(instruction.show(), 2)) for inst_num, instruction in enumerate(self.instructions)])

    def flatten_sems(self):
        """   Presents a string of non-empty sem properties so that all potentially synonymous sigs have the same flatten_sems, regardless of the order of the STRUCs.

        Used for development groping only, not used in production
        """
        rep = ''
        for instruction in self.instructions:
            props = [label for label in instruction.substantive_properties if instruction.get_property_values(label)]
            if props:
                rep += 'INSTRUCTION:' + ' '.join(props) + '|'
            for schedule in instruction.schedules:
                props = [label for label in schedule.substantive_properties if schedule.get_property_values(label)]
                if props:
                    rep += 'SCHEDULE:' + ' '.join(props) + '|'
                for event in schedule.events:
                    props = [label for label in event.substantive_properties if event.get_property_values(label)]
                    if props:
                        rep += 'EVENT:' + ' '.join(props) + '|'
        if len(rep) > 1 and rep[-1] == '|':
            rep = rep[:-1]
        if len(rep) > 2:
            rep = '||' + rep + '||'
        return rep





    def pprint_strucs(self):
        result = ''.join([struc.label if struc.is_unstructured_text() else '[' + struc.label + '|' + struc.quick_print_struc() + ']'  for struc in self.strucs])
        return result


    def __str__(self):
        non_blank_strucs = [struc for struc in self.strucs if not struc.label.isspace()]
        return '\n' + '\n'.join(['--Struc #%d: \n%s' % (struc_num, padlines(struc.__str__(), 5)) for struc_num, struc in enumerate(non_blank_strucs)])

    def show(self):
        row_1 = 'Raw Sig:       %s' % self.sig.raw_sig
        row_2 = 'Labels:        %s' % self.flatten()
        row_3 = 'Current State: %s' % self.__str__()
        return row_1 + '\n' + row_2 + '\n' + row_3

class Sem(object):
    """ A shell for Instruction, Schedule, and AdminEvent objects.
    """

    substantive_properties = []             # Substantive properties of the Sem. May differ from sem type to sem type.

    @classmethod
    def is_valid_property(cls, label):
        label = label.lower()
        return label in cls.substantive_properties

    def get_property_values(self, label):
        if not self.__class__.is_valid_property(label):
            return None
        else:
            return self.__dict__.get(label.lower())

    def add_property_value(self, label, struc):
        label = label.lower()

        if not self.__class__.is_valid_property(label):
            return None

        prop = self.__dict__[label]

        if type(prop) == list:
            prop.append(struc)
        else:
            self.__dict__[label] = struc

        struc.accounted_for_by_sem = self

    def remove_property(self, label, struc = None):

        label = label.lower()
        if not self.__class__.is_valid_property(label):
            return None

        prop = self.__dict__[label]

        if type(prop) == list and struc and struc in prop:
            prop.remove(struc)
        else:
            self.__dict__[label] = None

        if struc:
            struc.accounted_for_by_sem = None


    def show_props(self):
        props = [(label, self.get_property_values(label)) for label in self.substantive_properties]
        prop_reps = []
        for (label, prop) in props:
            if not prop:
                continue
            if type(prop) == list:
                val = '|'.join([a_prop.quick_print_struc() for a_prop in prop])
            else:
                val = prop.quick_print_struc()
            txt = (label + ': ').ljust(15) + val
            prop_reps.append(txt)

        return ('Props:\n' + padlines('\n'.join(prop_reps), 2)) if prop_reps else 'No Props'

    def get_prop_strucs(self):
        """ Returns the list of all strucs that are direct properties of the sem or are in the list of direct properties of a sem.

        E.g., if sem == Schedule and schedule.duration has value, adds that value (i.e. struc). If schedule.freq has 2 elements, includes them both.
        """

        props = [(label, self.get_property_values(label)) for label in self.substantive_properties]
        prop_strucs = []
        for (label, prop) in props:
            if not prop:
                continue
            if type(prop) == list:
                prop_strucs += prop
            else:
                prop_strucs.append(prop)
        return prop_strucs

    def get_recursive_componenets(self):
        """ Returns a list of all strucs that are properties of self or are properties of Schedules and Events that are in self.

            First come direct properties of the events in the first schedule, then properties of the first schedule,
            then properties of events of the second schedule, then properties of the second schedule,
            etc, and finally direct properties of the instruction.

        Is overriden by Instruction, Schedule, Event,
        """

        return []

    def __str__(self):
        return self.show()



class Instruction(Sem):
    """ Semantic object representing semantics of various types of complete thoughts - presumably Instructions.

    The Instruction is independent of the order of presentation or surface properties of the words in the sig. One should be able to mechanically generate
    a list of English sentences that represent the meaning of the Instruction.

    The typical Instruction (most sigs are just that one thought) is DrugAdmin: an instruction to administer some medication on a certain schedule.

    Other types of Instructions include directives to make an appointment, "Rinse out mouth after use" (for steroidal inhalers such as Beclomethasone for asthma),
    etc.

    Multiple Instructions of type DrugAdmin are possible in one parse if they refer to different drugs or different administration processes.

    parse               pointer to the parse where the Instruction is defined
    instruction_type    String: name of the type of instruction (typically, drugadmin)


    """

    substantive_properties = []             # Substantive properties of the Instruction. May differ from Instruction type to Instruction type.

    def __init__(self, parse, instruction_type):
        self.parse = parse
        self.instruction_type = instruction_type

    def get_coordinate_in_sem(self):
        """ Returns location of this instruction in the list of parse.instructions
        """
        coord_in_sem_list = self.parse.instructions.index(self)
        return (coord_in_sem_list, )

    def is_drugadmin_instruction(self):
        return self.instruction_type == 'drugadmin'

    def minimum_expected_number_of_atoms(self):
        """ Counts the number of events in the subparts of the Instruction. If not AdminEvent, returns 1.

            Used to estimate the minimal number of atoms that the Instruction could possibly be transduced to for quality of version calcs.
        """

        if not self.is_drugadmin_instruction():
            return 1

        count = 0
        count_from_substantive_properties = sum([1 for prop in self.substantive_properties if self.get_property_values(prop)])
        count += count_from_substantive_properties
        if not self.schedules:
            return count

        for sched_num, schedule in enumerate(self.schedules):
            if sched_num != 0:
                count += 1      # Each additional schedule usually requires at least 2 atoms "Then:" and the actual schedule substance.
            count_from_substantive_properties = sum([1 for prop in schedule.substantive_properties if schedule.get_property_values(prop)])
            count += max(count_from_substantive_properties - 1, 0)  # Expect at least 1 schedule property (e.g. Freq) to be covered by the first core sentence.
            if not schedule.events:
                count += 1      # If no events (and hence no core), add 1 to count_from_substantive_properties because we should expect 1 Atom for each schedule property
            else:
                count += (1 + 2 * (len(schedule.events) - 1))   # 1 for the first event and 2 for each subsequent event (because we need to add "Also:" before them).
        return count




    def show(self):
        return '\n' + '\n'.join(['--Prop %s: \n%s' % (label, padlines(self.get_property_values(label).__str__(), 5)) for label in self.substantive_properties if self.get_property_values(label)])

class DrugAdmin(Instruction):
    """ The main type of Instruction structure.

    schedules           List of Schedule objects. Typically used for tapering: each schedule specified the freq and duration of one tapering regime.

    indication          List: list of pointers to Indication(s) strucs.

    as_needed           Struc: pointer to AsNeeded() struc

    see_instructions    Struc: pointer to SeeInstructions() struc

    miscellaneous       Struc: pointer to Miscellaneous() struc

    primary_directive   Struc: pointer to an extracted or deduced for classification purposes from AdminEvents. Represents the directive ('take', 'inject', etc.) that dominates the meaning of the Sig.

    primary_form        Struc: pointer to extracted or deduced for classification purposes from AdminEvents. Represents the form ('tablet', 'applicatorful', etc.), if one is given or implied, for the main directive in the Sig.

    primary_route       Struc: pointer to extracted or deduced for classification purposes from AdminEvents.



    """

    substantive_properties = ['indication', 'as_needed', 'see_instructions',  'miscellaneous']

    def __init__(self, parse):
        Instruction.__init__(self, parse = parse, instruction_type = 'drugadmin')
        self.schedules = []
        self.indication = []
        self.as_needed = None
        self.see_instructions = None
        self.miscellaneous = None
        self.primary_directive = None
        self.primary_form = None
        self.primary_route = None


    def add_new_schedule(self):
        """ Adds a new empty schedule to the Instruction, and returns the new schedule  """

        schedule = Schedule(self)
        self.schedules.append(schedule)

        return schedule

    def get_recursive_componenets(self):
        """ Returns a list of all strucs that are properties of self or are properties of Schedules and Events that are in self.
            First come direct properties of the events in the first schedule, then properties of the first schedule,
            then properties of events of the second schedule, then properties of the second schedule,
            etc, and finally direct properties of the instruction.
        """

        components = []
        for schedule in self.schedules:
            components.extend(schedule.get_recursive_componenets())

        direct_components = self.get_prop_strucs()
        components.extend(direct_components)
        return components


    @staticmethod
    def is_valid_sem_property(label):
        """ Returns True if a struc with this label can be a property at any sem level: at Instruction(DrugAdmin), Schedule, or AdminEvent level.  """

        if not label:
            return False
        elif DrugAdmin.is_valid_property(label) or Schedule.is_valid_property(label) or AdminEvent.is_valid_property(label):
            return True
        else:
            return False


    def show(self):
        #props = ', '.join([label for label in self.substantive_properties if self.get_property_values(label)])
        props = self.show_props()
        schedules = '\n' + '\n'.join(['--Sched #%d: %s' % (sched_num, padlines(sched.show(), 0)) for sched_num, sched in enumerate(self.schedules)])
        schedules = padlines(schedules, 0) if self.schedules else ''
        header = 'Instruction type: %s with %d schedules' % (self.instruction_type, len(self.schedules) if 'schedules' in self.__dict__ else 0)
        return header + '\n' + padlines(props, 0) + schedules



class Schedule(Sem):
    """ A semantic structure for one schedule. E.g. "take 2 tabs in the AM and 1 in the PM", or

    instruction     pointer to the Instruction where the schedule belongs in the Instruction.schedules list

    events          List: List AdminEvent objects

    duration        Struc: duration struc

    stop_condition  Struc: e.g. till_gone

    maxdose         Struc: e.g. max 8 per 24 hrs.

    as_directed     Struc: AS_DIRECTED struc

    discard_remainder   Struc: DISCARD_REMAINDER

    repeat          Struc: REPEAT

    taper           Struc: TAPER

    freq            List of FREQ strucs

    periodicity     List of PERIODICITY strucs

    calendar_event  List: list of pointers to Calendar_Event descriptor Strucs.


    """

    substantive_properties = ['duration', 'stop_condition', 'maxdose', 'as_directed', 'discard_remainder', 'repeat', 'taper', 'freq', 'periodicity', 'calendar_event']


    def __init__(self, instruction):
        self.instruction = instruction
        self.events = []
        self.duration = None
        self.stop_condition = None
        self.maxdose = None
        self.as_directed = None
        self.discard_remainder = None
        self.repeat = None
        self.taper = None
        self.freq = []
        self.periodicity = []
        self.calendar_event = []

    def add_new_AdminEvent(self):
        """ Adds a new empty AdminEvent to the Schedule, and returns the new event    """

        event = AdminEvent(self)
        self.events.append(event)

        return event

    def get_coordinate_in_sem(self):
        """ Returns an ordered pair of integers representing the location of this schedule in the instruction structure.

        The first coordinate is the position of this schedule in the list schedules of the parent Instruction.
        The second coordinate is the location of the instruction on the parse.instructions list.
        """

        parent = self.instruction
        coord_in_schedules_list = parent.schedules.index(self)
        parent_coord_in_sems = parent.get_coordinate_in_sem()[0]
        return (coord_in_schedules_list, parent_coord_in_sems)

    def get_recursive_componenets(self):
        """ Returns a list of all strucs that are properties of self or are properties of Events that are in self.
            First come direct properties of the events in the schedule, then properties of the  schedule.
        """

        components = []
        for event in self.events:
            components.extend(event.get_recursive_componenets())

        direct_components = self.get_prop_strucs()
        components.extend(direct_components)
        return components

    def remove_event(self, event):
        """ Removes event from the schedule, and assigns all strucs that are part of the event to struc.accounted_for_by_sem = None

        It is not enough to just go through the component strucs of the Event. Other strucs (e.g. AND_CONJ) may be accounted_for_by_sem by this Event
        but are not components of the Event.
        """
        all_strucs = self.instruction.parse.strucs
        for struc in all_strucs:
            if struc.accounted_for_by_sem == event:
                struc.accounted_for_by_sem = None
        self.events.remove(event)
        event.schedule = None

    def show(self):
        props = self.show_props()
        events = '\n' + '\n'.join(['--Event #%d: %s' % (event_num, padlines(event.show(), 0)) for event_num, event in enumerate(self.events)])
        events = padlines(events, 0) if self.events else ''
        header = 'Sched with %d events' % (len(self.events))
        return header + '\n' + padlines(props + events, 2)

class AdminEvent(Sem):
    """ A semantic structure for one drug administration event. E.g. "1 tab by mouth for pain"

    directive       Pointer to Directive Struc
    dose            Pointer to Dose struc
    route           Pointer to Route
    site            Pointer to Site
    vehicle         Pointer to Vehicle (e.g. "via nebulizer")
    dir_attribute   List: list of pointers to DirAttribute (modifiers of Directive, such as "sparingly" or "take as 1 dose"
    anaphora        Pointer to Anaphora struc (e.g. "this medicine" -- used for Dictonary processing)
    timing          List: list of pointers to Timing descriptor Strucs.
    substrate       Pointer to Substrate
    specific_day    Pointer to Specific_Day

    """

    substantive_properties = ['directive', 'dose', 'route', 'site', 'vehicle', 'dir_attribute', 'anaphora', 'timing', 'substrate', 'specific_day']

    def __init__(self, schedule):
        self.schedule = schedule
        self.directive = None
        self.dose = None
        self.route = None
        self.site = None
        self.vehicle = None
        self.dir_attribute = []
        self.anaphora = None
        self.timing = []
        self.substrate = None
        self.specific_day = None


    def get_coordinate_in_sem(self):
        """ Returns an ordered triple (ev, sch, se) of integers representing the location of this Admin Event in the Instruction structure.

        The first coordinate is the position of this event in the parent schedule's list of events.
        The second coordinate is the position of the parent schedule in the list schedules of the grand-parent instruction.
        The third coordinate is the location of the grand-parent instruction on the parse.instructions list.
        """

        parent = self.schedule
        coord_in_events_list = parent.events.index(self)
        parent_coord_in_sems = parent.get_coordinate_in_sem()
        return (coord_in_events_list, parent_coord_in_sems[0], parent_coord_in_sems[1])

    def get_recursive_componenets(self):
        """ Returns a list of all strucs that are properties of self.    """

        components = self.get_prop_strucs()
        return components

    def show(self):
        return '\n' + padlines(self.show_props(),2)



class Struc(object):
    """
    A template for subcategorizing specific constituent structures.  All Strucs are listed in sequential order in which they are
    encountered in parse.strucs for each parse in sig.parses.

    label           is all-cap name for the type of constituent structure (e.g. SITE, NUM, DIR, etc.).
                    Unstructured (unparsed) chunks are labeled as themselves (i.e. lower case substrings of raw sig text).

    constituents    is a list of words/strucs this label expands to.

    value           Optional: Value is a specification of the instance of the label in question. E.g., if struc.label = 'FORM',
                    struc.value could be 'tablet', or if struc.label = 'NUM', struc.value could be 1.5, while
                    for a complex struc "one and a half tablets", the value is the list of values of
                    constituents:[1.5 tablet].

                    For unstructured text strucs, value is None.
    rules_used      List: For debugging: List of names of rules used to create/modify this structure

    accounted_for_by_sem        Pointer to a semantic structure (Instruction or Schedule or AdminEvent) that accounts for this information, if any.

    Two structures are equal if their labels and proper properties (those specified in proper_property_names) are point-wise equal.
    """

    def __init__(self, label, constituents = None):
        self.label = label
        self.constituents = constituents if constituents else []
        self.rules_used = []
        self.accounted_for_by_sem = None

    def is_valid_struc(self):
        """" Can be overriden with specific constraints. """
        return True

    def is_unstructured_text(self):
        """ A struc that is a substring of the original sig text, i.e. it has not been categorized.

            Since original sig text is lower-cased, while all labels for new strucs are upper case,
            a struc  is_unstructured_text() if it doesn't have any capital letters in its label.
        """

        return not self.label.isupper()

    def is_unparsed_non_punctuation_txt(self):
        """ The struc has not been parsed into a regular struc nor is it mere spaces or punctuation.        """

        pattern_non_label_text = re.compile('[a-z0-9]')
        if pattern_non_label_text.search(self.label):
            return True
        else:
            return False

    def is_space_or_punctuation_only(self):
        """ The struc is mere spaces or punctuation.

            Used in test_sems.py unit test.
        """

        pattern_space_or_punct = re.compile(r'[\s,\.:;\!\?\)\(\\\-]+')
        if pattern_space_or_punct.match(self.label):
            return True
        else:
            return False

    @property
    def property_names(self):
        props = sorted(list(set(self.__dict__.keys()) - set(['constituents', 'accounted_for_by_sem'])))
        return props

    @property
    def proper_property_names(self):
        props = sorted(list(set(self.__dict__.keys()) - set(['constituents', 'rules_used', 'accounted_for_by_sem'])))
        return props

    @property
    def substantive_properties(self):
        props = sorted(list(set(self.__dict__.keys()) - set(['constituents', 'rules_used', 'accounted_for_by_sem', 'label'])))
        return props

    def pointwise_equal(self, other):
        """ Two strucs are equal if they are pointwise equal (recursively for properties other than self.constituents and self.rules_used).

            We might want to loosen this requirement for some subclasses if we feel some properties are not required for identity.
        """

        if self.proper_property_names != other.proper_property_names:
            return False

        if self.label != other.label:
            return False

        for prop_name in self.proper_property_names:
            prop = self.__dict__[prop_name]
            other_prop = other.__dict__[prop_name]
            # all props are of type Struc or Int, String, or List (but in case of List it is never a list of Strucs)
            if isinstance(prop, Struc):
                if prop.pointwise_equal(other_prop):
                    pass
                else:
                    return False
            elif type(prop) == list:
                if set(prop) != set(other_prop):
                    return False
            else:
                if prop != other_prop:
                    return False

        return True



    def __str__(self):
        """
        The DOSE struc, for example ("2-4 tablets") will be represented as a string of
        <DOSE>
            quant:  <QUANT>
                        num_type:   range
                        low:        <QUANT>
                                        value: 2
                                        num_type: INT
                                    </QUANT>
                        high:       <QUANT>
                                        value: 4
                                        num_type: INT
                                    </QUANT>
                    </QUANT>
            form:   <FORM>
                        value: tablet
                        plurality: plural
                    </FORM>
        </DOSE>

        """

        if self.is_unstructured_text():
            return self.label
        else:
            result = ''
            tag = '<' + self.label + '>'
            closing_tag = '</' + self.label + '>'
            result += tag

            max_prop_name_len = max([len(prop_name) for prop_name in self.property_names] + [7])
            offset_len = max_prop_name_len + 7

            for prop_name in self.property_names:
                if prop_name == 'label':
                    continue
                prop = self.__dict__[prop_name]
                prop_repr = str(prop) if str(prop) != '' else 'None'
                prop_repr_lines = prop_repr.splitlines()
                first_line_prefix = ('\n' + ' ' * 4 + prop_name + ':' + ' ').ljust(offset_len)
                result += first_line_prefix + prop_repr_lines[0]
                for line in prop_repr_lines[1:]:
                    result += '\n' + ' ' * (offset_len - 1) + line
            if self.label == 'QUANT' and 'value' not in self.property_names:
                result += ('\n' + ' ' * 4 + 'value:' + ' ').ljust(offset_len) + str(self.value)
            if self.label == 'QUANT' and 'value_repr' not in self.property_names:
                result += ('\n' + ' ' * 4 + 'value_repr:' + ' ').ljust(offset_len) + self.value_repr
            if 'constituents' in self.__dict__:
                result += ('\n' + ' ' * 4 + 'const.:' + ' ').ljust(offset_len) + self.pretty_print_constituents()

            result += ('\n' + ' ' * 4 + 'key_set:' + ' ').ljust(offset_len) + str(get_key(self))

            result += '\n' + closing_tag

            return result

    def __repr__(self):
        return self.__str__()
        pass

    def show_accounted_for_by_sem(self, omit_coords = None):
        """ If this structure is accoutned for by a SEM, show which SEM and where (e.g. in the 3rd Schedule 1st event).

        Returns '|N|' if not accounted for by a SEM.
        Otherwise returns:
        if accounted for by an AdminEvent sem, returns '|AE_N_M_K|' where
            N = position of AdminEvent in sched.events,
            M = position of Schedule in instruction.schedules
            K = poistion of Instruction in parse.instructions
        if accounted for by an Schedule sem, returns '|SCH_M_K|' where
            M = position of Schedule in instruction.schedules
            K = poistion of Instruction in parse.instructions
        if accounted for by DrugAdmin Instruction, returns '|INSTDA_K|' where
            K = poistion of DrugAdmin Instruction in parse.instructions

        omit_coords         T/F. If True, only returns '|Y|' if struc is accounted for in an Instruction, omitting coordinates.
        """

        if self.is_unstructured_text():
            return ''

        sem_location = self.accounted_for_by_sem
        if not sem_location:
            return '|N|'
        elif omit_coords:
            return '|Y|'

        if isinstance(sem_location, AdminEvent):
            sem_type =  'AE'
        elif isinstance(sem_location, Schedule):
            sem_type =  'SCH'
        elif isinstance(sem_location, DrugAdmin):
            sem_type =  'INSTDA'

        coords = sem_location.get_coordinate_in_sem()

        return sem_type + '_' + '_'.join(str(coord) for coord in coords)

    def is_semantically_incompatible_with_given_sem(self, sem):
        """ Returns True if the struc can't be added to this sem as that sem as a property because that property is already taken (if only one slot per
            property is allowed, e.g. only one Route for AdminEvent.route) or, if the sem property is a list, if this struc contradicts
            one of the strucs assigned to that property on that sem.

            Returns True if the struc.label is not a valid property for that sem.

            This method gets overridden for specific Strucs for which the sem can have multiple compatible values, e.g. Timing and Dir_Attribute for AdminEvent or
            Indication for DrugAdmin or Freq for Schedule.
        """

        # Assuming that only unique property value is allowed for this sem property (otherwise should be overridden for list properties).
        if sem.get_property_values(self.label):
            return True
        else:
            # This struc can be added to this sem -- it is compatible.
            return False


    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        0   unacceptable
        1   not right, but acceptable if no other choices
        2   acceptable but defective
        3   acceptable but not perfect
        4   perfect match

        To be overridden by each Struc
        """

        raise Exception('Did not yet implement match_dictionary() for %s' % self.label)

    @staticmethod
    def match_dictionary_calc_multiple_scores(scores_list):
        """ Used in match_dictionary() calculations when there are multiple scores that need to be weighted to produce single match score.

            Allow some scores to be None, which means these are optional scores that are not applicable to the situation.
        """

        scores_list_pruned = [score for score in scores_list if score is not None]
        if not scores_list_pruned:
            return 4

        min_score = min(scores_list_pruned)
        average_score = sum(scores_list_pruned) / float(len(scores_list_pruned))

        if min_score == 0:
            return 0
        else:
            score = (min_score + average_score) / 2.0
        score = round(score, 2)
        return score


    def get_numerical_map_to_dict(self, dict_struc):
        """ For numerical variables in the dictionary such as <<NUM_1>> or <<DATE_0>> or <<TIME_0>>, returns a list of pairs
            [(var1, string repr of value1), (var2, value2), ...]. E.g., ('<<NUM_1>>', '1.5'), or ('<<DATE_0>>', '2/3/2012')

            For example, if struc == DOSE and DOSE.quant is a range 1 to 4, and in the dict DOSE is <<NUM_2>> to <<NUM_3>>,
            Then dose.get_numerical_map_to_dict() returns [('<<NUM_2>>', '1'), ('<<NUM_3>>', '4')]

            This is overridden for those Strucs that have a numerical component at some level, e.g. DOSE, DURATION, FREQ, Timing (e.g. "30 mins before breakfast") etc.
        """

        return []

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key

        For example, the key for a DOSE struc that represents "2-3 tablets" is DOSE/range-tablet. Maxdose of 12 tabs in 24 hours is MAXDOSE/tablet-hour.
        On the other hand, because Directive "Use" for dictionary atoms is compatible with any other directive, the key is empty when the directive is "use".
        The label is added by get_key_set() function because we only add the label at the top Sem level.


        To be overridden by each Struc
        """

        return None

    def pretty_print_constituents(self):
        """ Used for string representation of Strucs and Parses """

        if 'constituents' in self.__dict__:
            return ' '.join([item if type(item) in (str, unicode) else item.pretty_print_constituents() for item in self.constituents if item is not None])
        else:
            return ''


    def quick_print_struc(self):
        """ Used for string representation of Strucs in the Instruction.

            To be overriden by individual Strucs
        """


        if self.is_unstructured_text():
            return self.label

        raise Exception('Didnt expect to do quickprint for %s' % self.label)

        val = '-'.join([str(self.__dict__[prop]) if type(self.__dict__[prop]) in (str, bool, int, float) else prop.__dict__[prop].quick_print_struc() for prop in self.substantive_properties])
        return val

    def pprint(self):
        """ Returns stirngs such as 'DOSE/2-tablets'      """

        return self.label + ('' if self.is_unstructured_text() else '/'  + self.quick_print_struc())



    class StrucConstraintError(Exception):
        """ Violation of a constraint on a specific Struc.

        E.g. trying to create an instance of a Struc with the value of a slot being out of range.
        """

        pass


class Form(Struc):
    """Pill-type, e.g. Tablet, tablespoonful, pack.

    label           "FORM"
    value           String, the value of pill-type, eg tablet
    plurality       String, singular (tablet), plural (tablets), or either (e.g. "tablet(s)")
    permissible_values              List of possible values. None is permitted, too
    plurality_permissible_values    Set {'singular', 'plural', 'plurality_either'}
    value_to_likely_directive       Dict that maps forms to typical directives using them. It is used to deduce the directive when the form is given and directive is not.
    """

    permissible_values = ["application", "applicatorful", 'bottle', 'can', 'capful', "capsule", "cc", "drop", "dropperful", "gram", "inhalation", "liter", "lozenge", "mcg",
                          "mg", "ml", "ounce", "pack", "packet", "packaging", "pad", "patch", "puff", "ring", "scoop", "spray", "suppository", "syringe",
                          "tablespoon", "tablet", "teaspoon", "unit", "vial"]

    plurality_permissible_values = set(('singular', 'plural', 'plurality_either'))


    value_to_likely_directive = {'application': 'apply',
                                 'applicatorful': 'insert',
                                 'bottle': 'drink',
                                 'can': 'drink',
                                 'capful': 'mix',
                                 'capsule': 'take',
                                'cc': 'inject',
                                'drop': 'instill',
                                'dropperful': 'give',
                                'gram': 'apply',
                                'inhalation': 'use',
                                'lozenge': 'dissolve',
                                'mcg': 'inject',
                                'mg': 'inject',
                                'ml': 'take',
                                'ounce': 'rinse',
                                'packaging': 'unwrap',
                                'pack': 'use',
                                'packet': 'mix',
                                'patch': 'apply',
                                'puff': 'inhale',
                                'ring': 'insert',
                                'scoop': 'mix',
                                'spray': 'use',
                                'suppository': 'insert',
                                'syringe': 'inject',
                                'tablespoon': 'take',
                                'tablet': 'take',
                                'teaspoon': 'take',
                                'unit': 'inject',
                                'vial': 'inhale'}

    def __init__(self, form_name, plurality, constituents):
        Struc.__init__(self, label = 'FORM', constituents = constituents)
        self.value = form_name
        self.plurality = plurality if plurality else 'plurality_either'

        if debug and not self.is_valid_struc():
            error_msg = ('\n\nError specifying form value as -->%s<-- or of plurality of Form as -->%s<--. \nPermissible values of form_name are %s\nPermissible values of plurality of form are %s' %
                                       (form_name, plurality, ', '.join(self.permissible_values), ', '.join(list(self.plurality_permissible_values))))
            raise self.StrucConstraintError(error_msg)


    def is_valid_struc(self):
        value_valid = self.value is None or self.value in self.permissible_values
        plurality_valid = self.plurality in self.plurality_permissible_values
        return value_valid and plurality_valid

    def pointwise_equal(self, other):
        """ Two froms are equal if they have the same value, but not necessarily plurality (because it is often omitted or grammaticall incorrect)
        """
        return self.value == other.value

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.
        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        if self.value == 'unit':
            value = ''
        elif self.value in ('spray', 'puff'):
            value = 'spray|puff'
        else:
            value = self.value
        return value

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        0   unacceptable
        1   not right, but acceptable if no other choices
        2   acceptable but defective
        3   acceptable but not perfect
        4   perfect match
        """

        if self.value == dict_struc.value:
            value_score = 4
        elif self.value in ('puff', 'spray') and dict_struc.value in ('puff', 'spray'):
            value_score = 3
        elif dict_struc.value == 'unit' and self.value not in ('cc', 'gram', 'mcg', 'mg', 'liter', 'ounce'):
            value_score = 1
        else:
            return 0

        dict_rules = ', '.join(dict_struc.rules_used)
        self_rules = ', '.join(self.constituents)
        if 'deduced' in dict_rules and 'deduced' not in self_rules:
            value_score = min(value_score, 2)

        if self.plurality == dict_struc.plurality:
            plurality_score = 4
        elif dict_struc.plurality == 'plurality_either':
            plurality_score = 3.75
        else:
            plurality_score = 2

        scores = [value_score, plurality_score]
        return Struc.match_dictionary_calc_multiple_scores(scores)


    def quick_print_struc(self):
        return str(self.value) + ('s' if self.plurality != 'singular' else '')

class Vehicle(Struc):
    """ Vehicle of administration, e.g. "via nebulizer"

    value               String, e.g. nebulizer, inhaler
    preposition         String, e.g. via (nebulizer), in (inhaler), with (inhaler)
    permissible_values  List
    value_to_likely_directive   Dictionary mapping vehicle.value to possible directive values (as txt).
    """

    permissible_values = ('nebulizer', 'inhaler')
    value_to_likely_directive = {'nebulizer': 'inhale'}


    def __init__(self, value, preposition, constituents):
        Struc.__init__(self, label = 'VEHICLE', constituents = constituents)
        self.value = value
        self.preposition = preposition


    def is_valid_struc(self):
        test_result = self.value in self.permissible_values
        return test_result

    def pointwise_equal(self, other):
        """ Two vehicles are equal if they have the same value, but not necessarily preposition (because it is often omitted or grammaticall incorrect)
        """
        return self.value == other.value

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        """

        if self.value == dict_struc.value:
            return 4
        else:
            return 0

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        return self.value

    def quick_print_struc(self):
        return ((self.preposition + '-') if self.preposition else '') + self.value

class Quant(Struc):
    """ A number or a range of numbers.

    label           QUANT
    num_type        - one of ord|int|frac|decimal|range|var
                        where:
                            frac            The floating number is given in the sig as a fraction (e.g. "one and one half)".
                                            We use a subclass FracQuant to represent it then to keep the info about numerator and
                                            denominator, and to be able to provide alternative reprsentations as string depending whether
                                            it is part of a range or self-standing (e.g. "Take 1/2 (0.5) tabs" vs. "Take 1/2 - 2 tabs" (not appropriate to have parenthetical here).

                            decimal         The float is given in the sig as decimal "1.5" or "3.75ml"
                            ord             an ordinal, e.g first, second, 3rd, fifth, etc.
                            range           a range of numbers (int or float or ord), or a disjunction of numbers ("2 or 3 tabs").
                                            We use RangeQuant subclass to represent it.
                            var             Variable -- a special designation for dictionary variables such as <<NUM_1>> or <<DATE_0>> or <<TIME_0>>

    value           - If not range or var, then numerical value (int or float type), e.g. float(1)/3 or 1.5 or 400
                      If range, then the value is "range". If var, the value is 'var'
    value_repr      - a string with best representation of the value in numerical terms, e.g. "1+1/3" or "1+1/2" (1.5)"

    """

    def __init__(self, num_type, value = None):

        Struc.__init__(self, label = 'QUANT')
        self.num_type = num_type
        if value is not None:
            self.value = value
            self.value_repr = str(value)

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        0   unacceptable
        1   not right, but acceptable if no other choices
        2   acceptable but defective
        3   acceptable but not perfect
        4   perfect match
        """

        if self.num_type == dict_struc.num_type == 'range':
            low_score = self.low.match_dictionary(dict_struc.low)
            high_score = self.high.match_dictionary(dict_struc.high)
            if self.range_type == dict_struc.range_type:
                range_score = 4
            elif dict_struc.range_type == 'to':
                range_score = 3
            else:   # self is range of type "TO" (e.g. 2 to 4 pills) while dict_struc is of range type "OR" (e.g "1 or 2 pills").
                    # this will only work if self.low = 1 and self.high = 2 or self.low = 1/2 and self.high = 1
                if self.low.value == 1 and self.high.value == 2 and (dict_struc.low.num_type == 'var' or (dict_struc.low.value == 1 and dict_struc.high.value == 2)):
                    range_score = 3
                elif self.low.value == 0.5 and self.high.value == 1 and (dict_struc.low.num_type == 'var' or (dict_struc.low.value == 0.5 and dict_struc.high.value == 1)):
                    range_score = 3
                else:
                    range_score = 0
            return min(low_score, high_score, range_score)
        elif self.num_type != 'range' and dict_struc.num_type == 'var':
            if dict_struc.var_type == 'num':
                return 4
            else:
                return 0
        # The dictionary sometimes has fixed values for numbers, not variables. These are usually 1/2, 1, and 2
        elif self.value == dict_struc.value:
            return 4
        else:
            return 0

    def get_numerical_map_to_dict(self, dict_struc, omit_decimal_in_parenthesis = None):
        """ For numerical variables in the dictionary such as <<NUM_1>> or <<DATE_0>> or <<TIME_0>>, returns a list of pairs
            [(var1, string repr of value1), (var2, value2), ...]. E.g., ('<<NUM_1>>', '1.5'), or ('<<DATE_0>>', '2/3/2012')
        """

        if self.num_type == dict_struc.num_type == 'range':
            low = self.low.get_numerical_map_to_dict(dict_struc.low, omit_decimal_in_parenthesis = True)
            high = self.high.get_numerical_map_to_dict(dict_struc.high, omit_decimal_in_parenthesis = True)
            return low + high
        elif dict_struc.label == 'QUANT' and dict_struc.num_type == 'var' and dict_struc.var_type == 'num':
            var_num = dict_struc.var_number
            variable = '<<NUM_' + str(var_num) + '>>'
            if omit_decimal_in_parenthesis and self.num_type == 'frac':
                value = self.fraction_to_string(use_decimal_in_parenthesis = False)
            else:
                value = self.value_repr
            return [(variable, value)]
        elif dict_struc.label == 'QUANT' and dict_struc.num_type == 'var' and dict_struc.var_type == 'time':
            var_num = dict_struc.var_number
            variable = '<<TIME_' + str(var_num) + '>>'
            value = self.value_repr
            return [(variable, value)]
        else:
            return []

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Because Quant is not a Sem struc and doesn't exist in the atoms by itself, we only return here the value that is used in other Strucs
        containing Quant
        """

        if self.num_type == 'range':
            result = 'range'
        else:
            result = 'num'
        return result

    def quick_print_struc(self):
        return self.value_repr


class FracQuant(Quant):
    """
    For a frac "two and three-quarters" these properties take these values:

    integer         2
    numerator       3
    denominator     4
    value           2.75 (float)
    value_repr      "2+3/4 (2.75)" if use_decimal_in_parenthesis param to self.fraction_to_string() is True
    """

    def __init__(self, integer, numerator, denominator):
        self.integer = integer
        self.numerator = numerator
        self.denominator = denominator
        Quant.__init__(self, num_type = 'frac')


    def is_valid_struc(self, i = None, n = None, d = None):
        """ Checks validity of assumption that i = integer part of the fraction, n = numerator, d = denom. """

        if i is None:
            i = self.integer
        if n is None:
            n = self.numerator
        if d is None:
            d = self.denominator

        if not (type(i) == int and type(n) == int and type(d) == int):
            return False
        if not (i >= 0 and n > 0 and d > 0):     return False
        if n >= d: return False

        return True

    @property
    def value(self):
        if self.is_valid_struc(self.integer, self.numerator, self.denominator):
            return self.integer + float(self.numerator)/self.denominator
        elif debug:
            raise self.StrucConstraintError('Error creating Fraction with integer = %s, numerator = %s, denom = %s' %
                                       (str(self.integer), str(self.numerator), str(self.denominator)))
        else:
            return None

    @property
    def value_repr(self):
        return self.fraction_to_string(use_decimal_in_parenthesis = True)

    def fraction_to_string(self, use_decimal_in_parenthesis = None):
        """ Represents the number as string, e.g. "one and 2/3" -> "1+2/3", "31/2" -> "3+1/2"
            If the denominator is small and divides 100, we also express the answer as a decimal in
            parenthesis: "1 1/2" -> "1+1/2 (1.5)"
        """

        if self.integer == 0:
            i = ''
        else:
            i = str(self.integer) + '+'

        r = i + str(self.numerator) + '/' + str(self.denominator)

        if use_decimal_in_parenthesis and 100 % self.denominator == 0:       # i.e. the fraction is convertible to a small decimal that is best expressed as decimal
            r += ' (' + str(self.value) + ')'
        return r



class RangeQuant(Quant):
    """ A range of numbers (int or frac or decimal or ord) from low to high,  or a disjunction of numbers low or high ("2 or 3 tabs").
        Low and high are Quant strucs. If low or high is ord, then the other is also ord.

        range_type      String, to|or (e.g. "2 - 3 mg" is "to", "2 or 3 tabs" is "or"), a string
        low             Quant, lower bound
        high            Quant, upper bound
        value           'range'
        value_repr      str e.g. "2 to 3.5"
    """


    def __init__(self, range_type, low, high):

        Quant.__init__(self, num_type = 'range')
        self.low = low
        self.high = high
        self.range_type = range_type
        self.value = 'range'

        if debug and not self.is_valid_struc():
            raise self.StrucConstraintError('Error with range of low = %s and high = %s' %
                                       (str(low.value), str(high.value)))

    @property
    def value_repr(self):
        low = self.low
        high = self.high

        if low.num_type == 'frac':
            value_repr_low = low.fraction_to_string(use_decimal_in_parenthesis = False)
        else:
            value_repr_low = low.value_repr

        if high.num_type == 'frac':
            value_repr_high = high.fraction_to_string(use_decimal_in_parenthesis = False)
        else:
            value_repr_high = high.value_repr

        return value_repr_low + ' ' + self.range_type + ' ' + value_repr_high

    def is_valid_struc(self):
        try:
            if self.low.num_type == 'var' or self.high.num_type == 'var':
                test = self.low.num_type == self.high.num_type and self.low.var_type == self.high.var_type and self.low.var_number + 1 == self.high.var_number
                return test
            not_ranges = self.low.num_type != 'range' and self.high.num_type != 'range'
            # low and high can be either both ORD or neither of them is ORD. They can be
            # a free mix of INT, FRAC, or DECIMAL
            none_or_both_ord = not ( (self.low.num_type == 'ord') ^ (self.high.num_type == 'ord'))
            low_less_than_high = self.low.value < self.high.value
            return not_ranges and none_or_both_ord and low_less_than_high
        except:
            return False


class VarQuant(Quant):
    """ Variable -- a special designation for dictionary variables such as <<NUM_1>> or <<DATE_0>> or <<TIME_0>>

        var_type        String: num | date | time
        var_number      Integer, representing the sequence number of this variable in the dictionary string. E.g. for <<NUM_2>> var_number is 2.
        value           'var'
        value_repr      str e.g. "NUM_1", "TIME_0"
    """


    def __init__(self, var_type, var_number):

        Quant.__init__(self, num_type = 'var')
        self.var_type = var_type
        self.var_number = var_number
        self.value = 'var'


    @property
    def value_repr(self):
        return self.var_type.upper() + '_' + str(self.var_number)


class Date(Struc):
    """ The Struc of date type. Initialized with strings representing potential day, month, year. Represented as yyyy/mm/dd.

        day         Could be None, in which case the date is just yyyy/mm. O/w padded to 2 digits.
        month       Padded to 2 digits
        year        Corrected to 4 digits if needed.
        value       Represented as "yyyy/mm/dd"
    """


    def __init__(self, day, month, year):

        Struc.__init__(self, label = 'DATE')
        self.day = day
        self.month = month
        self.year = year

    @property
    def value(self):
        if not self.day:
            day = None
        else:
            day = ('0' + str(self.day))[-2:]
        month = ('0' + str(self.month))[-2:]
        year = ('20' + str(self.year))[-4:]
        if not day:
            value = year + '/' + month
        else:
            value = year + '/' + month + '/' + day

        return value

    def is_valid_struc(self):

        try:
            if self.day:
                day = int(self.day)
            else:
                day = None
            month = int(self.month)
            year = int(self.year)
        except ValueError:
            return False

        current_year = datetime.datetime.now().year
        if year < 100:
            year += 2000
            self.year = str(year)

        if month > 12 or month < 1 or year < 1990 or year >  current_year + 10:
            return False
        if type(day) == int and (day < 1 or day > 31):
            return False
        else:
            return True


    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        """

        if dict_struc.label == 'QUANT' and dict_struc.num_type == 'var' and dict_struc.var_type == 'date':
            return 4
        else:
            return 0

    def get_numerical_map_to_dict(self, dict_struc):
        """ For numerical variables in the dictionary such as <<NUM_1>> or <<DATE_0>> or <<TIME_0>>, returns a list of pairs
            [(var1, string repr of value1), (var2, value2), ...]. E.g., ('<<NUM_1>>', '1.5'), or ('<<DATE_0>>', '2/3/2012')
        """

        if dict_struc.label == 'QUANT' and dict_struc.num_type == 'var' and dict_struc.var_type == 'date':
            var_num = dict_struc.var_number
            variable = '<<DATE_' + str(var_num) + '>>'
            value = self.value
            return [(variable, value)]
        else:
            return []



    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """
        return True

    def quick_print_struc(self):
        return self.value

class BloodPressureRatio(Struc):
    """ The Struc of type n/m which means blood pressure systolic/diastolic'

        systolic    Int: The upper number in the "frac"
        diastolic   Int: The lower number in the "frac"
        value       Str: "sys/dias"
    """

    def __init__(self, systolic, diastolic):
        """systolic and diastolic are integers"""

        Struc.__init__(self, label = 'BLOODPRESSURERATIO')
        self.systolic = int(systolic)
        self.diastolic = int(diastolic)

    @property
    def value(self):
        return str(self.systolic) + '/' + str(self.diastolic)

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        """
        return 4

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """
        return True

    def quick_print_struc(self):
        return self.value

class Directive(Struc):
    """ The verb

    value                           String: the verb
    permissible_values              List of permissible Directives.
    """

    permissible_values = ['apply', 'chew and swallow', 'chew', 'dissolve', 'drink', 'give', 'inhale', 'inject',
                          'insert', 'instill', 'mix', 'place', 'put', 'remove', 'rinse', 'shampoo', 'spray', 'stop', 'swallow',
                          'swish', 'swish and spit', 'swish and swallow', 'take', 'test', 'unwrap and insert', 'unwrap', 'use']

    directive_to_likely_form = {'apply': 'patch',
                                'chew and swallow': 'tablet',
                                'chew': 'tablet',
                                'dissolve': 'tablet',
                                'give': 'tablet',
                                'inhale': 'puff',
                                'inject': 'unit',
                                'insert': 'applicatorful',
                                'instill': 'drop',
                                'mix': 'gram',
                                'place': 'drop',
                                'remove': 'patch',
                                'spray': 'puff',
                                'take': 'tablet',
                                'unwrap': 'packaging',
                                'unwrap and insert': 'suppository',
                                'use': ''}

    def __init__(self, value):
        Struc.__init__(self, label = 'DIRECTIVE')
        self.value = value

    def is_valid_struc(self):
        test_result = self.value in  self.permissible_values
        return test_result


    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        0   unacceptable
        1   not right, but acceptable if no other choices
        2   acceptable but defective
        3   acceptable but not perfect
        4   perfect match
        """

        if self.value == dict_struc.value:
            score = 4
        elif self.value == 'give' and dict_struc.value == 'take':
            score = 2
        elif dict_struc.value == 'use':
            score = 2
        else:
            score = 0

        return score

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.
            Because "use" is compatible with everything, the key for "use" is empty.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        if self.value == 'use':
            value = ''
        elif self.value == 'give':
            value = 'take'
        else:
            value = self.value

        return value

    def quick_print_struc(self):
        return self.value

class DirAttribute(Struc):
    """ Modifiers of the Directive, e.g. "take sparingly", "take as 1 dose", "test blood sugar", "per sliding scale"

    There can be multiple dir_attributes in a Sig (e.g. "test blood sugar per sliding scale"). To make sure that string_key() set
    of a sig matches all possible dictionary atoms, we only allow dir_attribute to have one string value, but group multiple
    dir_attributes in a list at the Admin_Event sem level.

    value                   Strings
    permissible_values      Set of strings of permissible values for self.values

    """

    permissible_values = set(['as_single_dose', 'blood_sugar', 'first_dose', 'liberally', 'lightly', 'once_only', 'per_sliding_scale', 'small_amount', 'sparingly', 'thinly', 'together'])


    def __init__(self, value):
        Struc.__init__(self, label = 'DIR_ATTRIBUTE')
        self.value = value

    def is_valid_struc(self):
        test_result = self.value in self.permissible_values
        return test_result


    def is_semantically_incompatible_with_given_sem(self, admin_event):
        """ Returns True if the dir_attribute structure is incompatible with at least one of the dir_attribute structures on the AdminEvent.dir_attribute list
            (which is a list of dir_attribute strucs). I.e. Self is incompatible with the dir_attribute list of an AdminEvent if Self can't be added to that AdminEvent.dir_attribute -- it signals
            a new AdminEvent.

        For now, this is a dummy returning False.
        """

        return False

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        """

        if self.value == dict_struc.value:
            score = 4
        else:
            score = 0

        return score

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        return self.value

    def quick_print_struc(self):
        return self.value

class Dose(Struc):
    """ Represents the quantity of form to take, e.g. "2-3 puffs"
    quant           Quant struc, e.g. a range "2 to 3"
    form            Form struc, e.g. tablet.
    up_to           Optional: True/False flag to signify if the Quant is meant to be upper limit. E.g. "upto 3 puffs"
    constituents    [quant, form]
    """

    def __init__(self, quant, form):
        Struc.__init__(self, label = 'DOSE', constituents = [quant, form])
        self.quant = quant
        self.form = form
        self.up_to = None


    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        0   unacceptable
        1   not right, but acceptable if no other choices
        2   acceptable but defective
        3   acceptable but not perfect
        4   perfect match
        """

        quant_score = self.quant.match_dictionary(dict_struc.quant)
        # Because for all quants in self besides ranges the quant_score is going to be 4 or 0, it unfairly dilutes the contribution of form to quality.
        # This is especially bothersome for choosing between things like "puff" or "spray" vs. "unit".
        if quant_score == 0:
            return 0
        elif self.quant.num_type in ('int', 'frac', 'decimal') and quant_score == 4:
            quant_score = None

        form_score = self.form.match_dictionary(dict_struc.form)
        if self.up_to and dict_struc.up_to:
            up_to_score = 4
        elif self.up_to:
            up_to_score = 1
        elif dict_struc.up_to:
            # self.up_to is None
            return 0
        else:
            # both are None
            up_to_score = None

        scores = [quant_score, form_score, up_to_score]
        return Struc.match_dictionary_calc_multiple_scores(scores)

    def get_numerical_map_to_dict(self, dict_struc):
        """ For numerical variables in the dictionary such as <<NUM_1>> or <<DATE_0>> or <<TIME_0>>, returns a list of pairs
            [(var1, string repr of value1), (var2, value2), ...]. E.g., ('<<NUM_1>>', '1.5'), or ('<<DATE_0>>', '2/3/2012')
        """

        found = self.quant.get_numerical_map_to_dict(dict_struc.quant)
        return found

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        if self.quant.num_type == 'range':
            quant = 'range-'
        else:
            quant = ''

        if self.form and self.form.value:
            form = self.form.string_key()
        else:
            form = ''


        return quant + form

    def quick_print_struc(self):
        up_to = 'upto_' if self.up_to else ''
        return up_to + self.quant.quick_print_struc() + ('-' + self.form.quick_print_struc() if self.form else '')


class Substrate(Struc):
    """ Represents the substrate that is used to mix or to drink the dose with.

        E.g.
            mix 17 gms with 8 oz of water
            take 1 tablet twice daily with a glass of water

    value           String: water | liguid
    dose            Optional: Dose struc, e.g. 8 oz
    """

    permissible_values = ('liquid', 'water', 'gatorade', 'diluent')


    def __init__(self, value, dose = None):
        Struc.__init__(self, label = 'SUBSTRATE')
        self.value = value
        self.dose = dose

    def is_valid_struc(self):
        test_result = (self.value in self.permissible_values)
        return test_result


    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        0   unacceptable
        1   not right, but acceptable if no other choices
        2   acceptable but defective
        3   acceptable but not perfect
        4   perfect match
        """

        if self.value == dict_struc.value:
            value_score = 4
        elif dict_struc.value == 'liquid':
            value_score = 1
        else:
            return 0

        if self.dose and dict_struc.dose:
            dose_score = self.dose.match_dictionary(dict_struc.dose)
        elif self.dose is None and dict_struc.dose is None:
            dose_score = 4
        elif dict_struc.dose is None:
            dose_score = 2
        else:
            return 0

        scores = [value_score, dose_score]
        return Struc.match_dictionary_calc_multiple_scores(scores)


    def get_numerical_map_to_dict(self, dict_struc):
        """ For numerical variables in the dictionary such as <<NUM_1>> or <<DATE_0>> or <<TIME_0>>, returns a list of pairs
            [(var1, string repr of value1), (var2, value2), ...]. E.g., ('<<NUM_1>>', '1.5'), or ('<<DATE_0>>', '2/3/2012')
        """

        if self.dose and dict_struc.dose:
            dose = self.dose.get_numerical_map_to_dict(dict_struc.dose)
        else:
            dose = []

        return dose

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key() procedure.
        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        return None

    def quick_print_struc(self):
        if self.dose:
            dose = self.dose.quick_print_struc()
        else:
            dose = ''
        return self.value + ('/' + dose if dose else '')



class TimeUnit(Struc):
    """ Minute, Hour, Day, Week, Month

    label       "TIMEUNIT"
    value       String, second|minute|hour|day|week|month
    plurality   singular (hour), plural (hours)
    plurality_permissible_values    set {'singular', 'plural', 'plurality_either'}
    """

    permissible_values = ('second', 'minute', 'hour', 'day', 'week', 'month')
    plurality_permissible_values = set(('singular', 'plural', 'plurality_either'))

    def __init__(self, value, constituents, plurality = None):
        Struc.__init__(self, label = 'TIMEUNIT', constituents = constituents)
        self.value = value
        if not plurality:
            plurality = 'singular'
        self.plurality = plurality      #   values for plurality can be 'singular', 'plural', 'plurality_either'

        if debug and not self.is_valid_struc():
            error_msg = ('Error specifying plurality of TimeUnit as %s. Permissible values are %s' %
                                       (plurality, ', '.join(list(self.plurality_permissible_values))))
            raise self.StrucConstraintError(error_msg)


    def is_valid_struc(self):
        test_result = self.value in self.permissible_values and self.plurality in self.plurality_permissible_values
        return test_result

    def pointwise_equal(self, other):
        """ Two time units are equal if they have the same value, but not necessarily plurality (because it is often omitted or grammaticall incorrect)
        """

        return self.value == other.value

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        return self.value

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        0   unacceptable
        1   not right, but acceptable if no other choices
        2   acceptable but defective
        3   acceptable but not perfect
        4   perfect match
        """

        if self.value != dict_struc.value:
            return 0

        if self.plurality == dict_struc.plurality:
            plurality_score = 4
        elif dict_struc.plurality == 'plurality_either':
            plurality_score = 3
        else:
            plurality_score = 2

        return plurality_score

    def quick_print_struc(self):
        return self.value + ('s' if self.plurality != 'singular' else '')

class TimeInterval(Struc):
    """ Represents the time interval, e.g. 2 hours, 3 weeks.

    quant           QUANT struc
    time_unit       TIMEUNIT struc
    constituents    [quant, time_unit]
    """

    def __init__(self, quant, time_unit):
        Struc.__init__(self, label = 'TIMEINTERVAL', constituents = [quant, time_unit])
        self.quant = quant
        self.time_unit = time_unit

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        """

        quant_score = self.quant.match_dictionary(dict_struc.quant)
        time_unit_score = self.time_unit.match_dictionary(dict_struc.time_unit)

        scores = [quant_score, time_unit_score]
        return Struc.match_dictionary_calc_multiple_scores(scores)

    def get_numerical_map_to_dict(self, dict_struc):
        """ For numerical variables in the dictionary such as <<NUM_1>> or <<DATE_0>> or <<TIME_0>>, returns a list of pairs
            [(var1, string repr of value1), (var2, value2), ...]. E.g., ('<<NUM_1>>', '1.5'), or ('<<DATE_0>>', '2/3/2012')
        """

        found = self.quant.get_numerical_map_to_dict(dict_struc.quant)
        return found

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        return self.quant.string_key() + '-' + self.time_unit.string_key()

    def quick_print_struc(self):
        return self.quant.quick_print_struc() + '-' + self.time_unit.quick_print_struc()

class MaxDose(Struc):
    """ Represents the maximum dose one can take in a time interval, e.g. "do not exceed 4 tabs / 12 hour".

    dose            DOSE struc, but dose.form could be None if unspecified
    time_interval   TimeInterval struc
    directive       Directive struc (e.g. "do not INSTILL more than ..."
    constituents    [dose, time_interval, directive]
    """

    def __init__(self, dose, time_interval, directive):
        Struc.__init__(self, label = 'MAXDOSE', constituents = [dose, time_interval, directive])
        self.dose = dose
        self.time_interval = time_interval
        self.directive = directive

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        0   unacceptable
        1   not right, but acceptable if no other choices
        2   acceptable but defective
        3   acceptable but not perfect
        4   perfect match
        """

        dose_score = self.dose.match_dictionary(dict_struc.dose)
        if dose_score == 0:
            quant_score = self.dose.quant.match_dictionary(dict_struc.dose.quant)
            form_score = self.dose.form.match_dictionary(dict_struc.dose.form)
            if form_score == 0 and dict_struc.dose.form and dict_struc.dose.form.value == 'unit':
                form_score = 2
            dose_score = min(quant_score, form_score)
        time_interval_score = self.time_interval.match_dictionary(dict_struc.time_interval)
        if self.directive and dict_struc.directive:
            directive_score = self.directive.match_dictionary(dict_struc.directive)
        elif dict_struc.directive:
            if dict_struc.directive.value == 'use':
                directive_score = 4
            elif dict_struc.directive.value == 'take':
                directive_score = 2
            else:
                directive_score = 0
        elif self.directive:
            if self.directive.value == 'use':
                directive_score = 4
            elif self.directive.value == 'take':
                directive_score = 3
            else:
                directive_score = 2
        else:           # neither side has a directive
            directive_score = 4

        scores = [dose_score, time_interval_score, directive_score]
        return Struc.match_dictionary_calc_multiple_scores(scores)

    def get_numerical_map_to_dict(self, dict_struc):
        """ For numerical variables in the dictionary such as <<NUM_1>> or <<DATE_0>> or <<TIME_0>>, returns a list of pairs
            [(var1, string repr of value1), (var2, value2), ...]. E.g., ('<<NUM_1>>', '1.5'), or ('<<DATE_0>>', '2/3/2012')
        """

        dose = self.dose.get_numerical_map_to_dict(dict_struc.dose)
        if self.time_interval:
            time_interval = self.time_interval.get_numerical_map_to_dict(dict_struc.time_interval)
        else:
            time_interval = []
        return dose + time_interval

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        if self.directive:
            directive = self.directive.string_key()
        else:
            directive = None
        directive = (directive + '-') if directive else ''
        dose = self.dose.string_key()
        time_interval = self.time_interval.string_key()
        time_interval = ('-' + time_interval) if time_interval else time_interval

        return time_interval

    def quick_print_struc(self):
        directive = (self.directive.value + '-') if self.directive else ''
        return 'MaxDose:' + directive + self.dose.quick_print_struc() + '/' + self.time_interval.quick_print_struc()

class Duration(Struc):
    """ Represents the duration concept, e.g. "take for 3 days", "take for the next 3 days", "for 3 more days"

    time_interval   TimeInterval struc

    on_off          String (optional). Is set to 'on' if the duration is part of the front end of expression "on for a period then off for a period"
                    "APPLY ONE PATCH LEAVE ON FOR 12 HOURS THEN OFF FOR 12 HOURS"
                    Is set to 'off' if the duration is part of the back end of expression "on for a period then off for a period"
                    Only is set when time_interval.time_unit = 'hour' because that is the only case when the existing Sem structures (which assume that Duration is on
                    inter-day level) don't apply appropriately. See rule_on_before_off(), rule_off_after_on(), and process_special_duration_cases().

    offset          String (optional). Means that the duration starts not immediately but:
                    if offset = 'next', duration starts after the current time period ends. E.g. "take for the next 3 days" or "take for 3 more days"
                    if offset = 'first', duration starts immediately, e.g. "for the first 3 days:"

    up_to           True/False (optional). Used with, e.g. "apply for up to 12 hours."


    constituents    [TimeInterval]
    """

    def __init__(self, time_interval):
        Struc.__init__(self, label = 'DURATION', constituents = [time_interval])
        self.time_interval = time_interval
        self.on_off = None
        self.offset = ''
        self.up_to = None

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        0   unacceptable
        1   not right, but acceptable if no other choices
        2   acceptable but defective
        3   acceptable but not perfect
        4   perfect match
        """

        time_interval_score = self.time_interval.match_dictionary(dict_struc.time_interval)

        if self.offset == dict_struc.offset:
            offset_score = 4
        else:
            return 0

        if self.on_off == dict_struc.on_off:
            on_off_score = 4
        else:
            return 0

        if self.up_to and dict_struc.up_to:
            up_to_score = 4
        elif self.up_to:
            up_to_score = 1
        elif dict_struc.up_to:
            # self.up_to is None
            return 0
        else:
            # both are None
            up_to_score = None

        scores = [time_interval_score, offset_score, on_off_score, up_to_score]
        return Struc.match_dictionary_calc_multiple_scores(scores)

    def get_numerical_map_to_dict(self, dict_struc):
        """ For numerical variables in the dictionary such as <<NUM_1>> or <<DATE_0>> or <<TIME_0>>, returns a list of pairs
            [(var1, string repr of value1), (var2, value2), ...]. E.g., ('<<NUM_1>>', '1.5'), or ('<<DATE_0>>', '2/3/2012')
        """

        if self.time_interval:
            time_interval = self.time_interval.get_numerical_map_to_dict(dict_struc.time_interval)
        else:
            time_interval = []

        return time_interval

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        time_interval = self.time_interval.string_key()
        offset = self.offset
        offset = ('-' + offset) if offset else ''

        return time_interval + offset


    def quick_print_struc(self):
        on_off_prefix = (self.on_off + '-') if self.on_off else ''
        offset = (self.offset + '-') if self.offset else ''
        return  on_off_prefix + 'for-' + offset + self.time_interval.quick_print_struc()

class Freq(Struc):
    """ Represents frequency, i.e. quant per timeunit, e.g "2 times daily"
    quant           QUANT struc, e.g. a range "2 to 3"
    time_unit       TIMEUNIT struc, e.g. week
    up_to           Optional: Flag True/False. If True, means "up to Quant" e.g. "up to 3 times a day"
    constituents    [quant, time_unit]
    """

    def __init__(self, quant, time_unit):
        Struc.__init__(self, label = 'FREQ', constituents = [quant, time_unit])
        self.quant = quant
        self.time_unit = time_unit
        self.up_to = None

    def is_semantically_incompatible_with_given_sem(self, schedule):
        """ Returns True if this Freq struc is incompatible with at least one of the Freq structures on the schedule.freq list. I.e., the struc can't be
            added to this schedule -- it signals that a new Schedule needs to be started

            A Schedule can have several Frequencies associated with it. E.g. "take 2 tabs 3 times a day twice weekly"  which has intra-day and intra-week freq.

            These Freqs can't be in the same schedule if they share the same time_unit but not the same quant.
            We allow pointwise identical Freqs because there are lots of stupid repetitions like "take daily 2 times a day"
        """

        quant = self.quant
        time_unit = self.time_unit

        for freq in schedule.freq:
            if time_unit.pointwise_equal(freq.time_unit) and not quant.pointwise_equal(freq.quant):
                return True
        return False

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        """

        quant_score = self.quant.match_dictionary(dict_struc.quant)
        # Because for all quants in self besides ranges the quant_score is going to be 4 or 0, it unfairly dilutes the contribution of other factors (time_unit, up_to) to quality.
        if quant_score == 0:
            return 0
        elif self.quant.num_type in ('int', 'frac', 'decimal') and quant_score == 4:
            quant_score = None

        time_unit_score = self.time_unit.match_dictionary(dict_struc.time_unit)

        if self.up_to and dict_struc.up_to:
            up_to_score = 4
        elif self.up_to:
            # allow "up to 3 times/day" to match to "3 times/day" as a desperation
            up_to_score = 1
        elif dict_struc.up_to:
            # self.up_to is None
            return 0
        else:
            # both are None
            up_to_score = None

        scores = [quant_score, time_unit_score, up_to_score]
        return Struc.match_dictionary_calc_multiple_scores(scores)

    def get_numerical_map_to_dict(self, dict_struc):
        """ For numerical variables in the dictionary such as <<NUM_1>> or <<DATE_0>> or <<TIME_0>>, returns a list of pairs
            [(var1, string repr of value1), (var2, value2), ...]. E.g., ('<<NUM_1>>', '1.5'), or ('<<DATE_0>>', '2/3/2012')
        """

        quant = self.quant.get_numerical_map_to_dict(dict_struc.quant)
        return quant

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        quant = self.quant.string_key()
        time_unit = self.time_unit.string_key()

        return quant + '-' + time_unit

    def quick_print_struc(self):
        up_to = 'upto_' if self.up_to else ''
        return up_to + self.quant.quick_print_struc() + '/' + self.time_unit.quick_print_struc()

class Periodicity(Struc):
    """ Represents periodicity, e.g. every n'th quant timeunit, e.g "every 4-8 hours", "every other day", "every 3rd week"
    I.e. how many time units to wait between dosing events.

    quant           QUANT struc representing how many time units to wait, e.g. a range "2 to 3"
    time_unit       TIMEUNIT struc, e.g. week
    up_to           Optional: True/False flag. E.g. "up to every 8 hours"
    constituents    [quant, time_unit]
    """

    def __init__(self, quant, time_unit):
        Struc.__init__(self, label = 'PERIODICITY', constituents = [quant, time_unit])
        self.quant = quant
        self.time_unit = time_unit
        self.up_to = None


    def is_semantically_incompatible_with_given_sem(self, schedule):
        """ Returns True if this Periodicity struc is incompatible with at least one of the Periodicity structures on the schedule.periodicity list. I.e.,
            the struc can't be added to this schedule. It signals that a new Schedule needs to be started.

            A Schedule can have several Periodicities associated with it. E.g. "take 2 tabs every 8 hours every day", or more typically,
            "1 tab daily every 8 hours" which means "take 1 tab every 8 hours. Take it every day."

            These Periodicities can't be in the same schedule if they share the same time_unit but not the same quant.
            We allow pointwise identical Periodicities because there are lots of stupid repetitions like "take daily 2 times a day"
        """

        quant = self.quant
        time_unit = self.time_unit

        for struc in schedule.periodicity:
            if time_unit.pointwise_equal(struc.time_unit) and not quant.pointwise_equal(struc.quant):
                return True
        return False

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        """

        quant_score = self.quant.match_dictionary(dict_struc.quant)
        # Because for all quants in self besides ranges the quant_score is going to be 4 or 0, it unfairly dilutes the contribution of other factors (time_unit, up_to) to quality.
        if quant_score == 0:
            return 0
        elif self.quant.num_type in ('int', 'frac', 'decimal') and quant_score == 4:
            quant_score = None

        time_unit_score = self.time_unit.match_dictionary(dict_struc.time_unit)

        if self.up_to and dict_struc.up_to:
            up_to_score = 4
        elif self.up_to:
            # allow "up to every 4 hours" to match to "every 4 hours" as a desperation
            up_to_score = 1
        elif dict_struc.up_to:
            # self.up_to is None
            return 0
        else:
            # both are None
            up_to_score = None

        scores = [quant_score, time_unit_score, up_to_score]
        return Struc.match_dictionary_calc_multiple_scores(scores)

    def get_numerical_map_to_dict(self, dict_struc):
        """ For numerical variables in the dictionary such as <<NUM_1>> or <<DATE_0>> or <<TIME_0>>, returns a list of pairs
            [(var1, string repr of value1), (var2, value2), ...]. E.g., ('<<NUM_1>>', '1.5'), or ('<<DATE_0>>', '2/3/2012')
        """

        quant = self.quant.get_numerical_map_to_dict(dict_struc.quant)
        return quant

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        quant = self.quant.string_key()
        time_unit = self.time_unit.string_key()

        return quant + '-' + time_unit

    def quick_print_struc(self):
        up_to = 'upto_' if self.up_to else ''
        return up_to + 'every-' + self.quant.quick_print_struc() + '-' + self.time_unit.quick_print_struc()


class Repeat(Struc):
    """ Represents Repeat instruction, e.g.:
            "may repeat every 5 minutes up to 3 times",
            "take one tablet by mouth at onset of headache may repeat in 2 hours"
            "take 1 tablet at onset of migraine. may repeat once after 2 hours. max 10 mg/day."
            "1 tablet once daily and repeat 1 tablet 1 week later"
            "apply 1 patch once weekly as directed for 3 weeks, leave off for 1 week then repeat cycle"

    periodicity     PERIODICITY struc: Optional, representing how many time units to wait, e.g. "every 5 minutes"
    offset          TimeInterval struc: Optional, start repeating after this interval (e.g. may repeat after 2 hours)
    reps            QUANT struc: Optional, repeat this many times
    max_reps        QUANT struc: Optional, up to how many times to repeat (can't be both reps and max_reps)
    condition       String: Optional, condition for repeating (e.g. "if not improvement"
    constituents    [found text]
    """

    condition_permissible_values = ('no_relief')

    def __init__(self, periodicity, offset, reps, max_reps, condition):
        Struc.__init__(self, label = 'REPEAT')
        self.periodicity = periodicity
        self.offset = offset
        self.reps = reps
        self.max_reps = max_reps
        self.condition = condition

    def is_valid_struc(self):
        condition_test = self.condition in self.condition_permissible_values if self.condition else True
        reps_test = not (self.reps and self.max_reps)
        return condition_test and reps_test

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        """

        if self.periodicity and dict_struc.periodicity:
            periodicity_score = self.periodicity.match_dictionary(dict_struc.periodicity)
        elif self.periodicity or dict_struc.periodicity:
            return None
        else:
            periodicity_score = 4

        if self.offset and dict_struc.offset:
            offset_score = self.offset.match_dictionary(dict_struc.offset)
        elif self.offset or dict_struc.offset:
            return None
        else:
            offset_score = 4

        if self.reps and dict_struc.reps:
            reps_score = self.reps.match_dictionary(dict_struc.reps)
        elif self.reps and dict_struc.max_reps:     # "may repeat" somewhat naturally implies "up to", so we can allow "may repeat 3 times" to match "may repeat up to 3 times"
            reps_score = min(2, self.reps.match_dictionary(dict_struc.max_reps))
        elif self.reps or dict_struc.reps:
            return None
        else:
            reps_score = 4

        if self.max_reps and dict_struc.max_reps:
            max_reps_score = self.max_reps.match_dictionary(dict_struc.max_reps)
        elif self.reps and dict_struc.max_reps:
            max_reps_score = reps_score
        elif self.max_reps or dict_struc.max_reps:
            return None
        else:
            max_reps_score = 4

        if self.condition and dict_struc.condition:
            condition_score = (self.condition == dict_struc.condition)
        elif self.condition or dict_struc.condition:
            return None
        else:
            condition_score = 4

        scores = [periodicity_score, offset_score, reps_score, max_reps_score, condition_score]
        return Struc.match_dictionary_calc_multiple_scores(scores)

    def get_numerical_map_to_dict(self, dict_struc):
        """ For numerical variables in the dictionary such as <<NUM_1>> or <<DATE_0>> or <<TIME_0>>, returns a list of pairs
            [(var1, string repr of value1), (var2, value2), ...]. E.g., ('<<NUM_1>>', '1.5'), or ('<<DATE_0>>', '2/3/2012')
        """

        if self.periodicity and dict_struc.periodicity:
            periodicity = self.periodicity.get_numerical_map_to_dict(dict_struc.periodicity)
        else:
            periodicity = []

        if self.offset and dict_struc.offset:
            offset = self.offset.get_numerical_map_to_dict(dict_struc.offset)
        else:
            offset = []

        if self.reps and dict_struc.reps:
            reps = self.reps.get_numerical_map_to_dict(dict_struc.reps)
        elif self.reps and dict_struc.max_reps:
            reps = self.reps.get_numerical_map_to_dict(dict_struc.max_reps)
        else:
            reps = []

        if self.max_reps and dict_struc.max_reps:
            max_reps = self.max_reps.get_numerical_map_to_dict(dict_struc.max_reps)
        else:
            max_reps = []

        return periodicity + offset + reps + max_reps



    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        periodicity = ('every-' + self.periodicity.string_key()) if self.periodicity else ''
        offset = ('|after-' + self.offset.string_key()) if self.offset else ''
        if self.reps:
            reps = '|reps-' + self.reps.string_key()
        elif self.max_reps:
            reps = '|reps-' + self.max_reps.string_key()
        else:
            reps = ''
        condition = ('|condition-' + self.condition) if self.condition else ''

        return periodicity + offset + reps + condition

    def quick_print_struc(self):
        periodicity = self.periodicity.quick_print_struc() if self.periodicity else ''
        offset = ('/after-' + self.offset.quick_print_struc()) if self.offset else ''
        reps = ('/reps-' + self.reps.quick_print_struc()) if self.reps else ''
        max_reps = ('/max_reps-' + self.max_reps.quick_print_struc()) if self.max_reps else ''
        condition = ('/condition-' + self.condition) if self.condition else ''
        return periodicity + offset + reps + max_reps + condition



class Taper(Struc):
    """ Represents Taper instruction, e.g.:
            "take one capsule three times daily may increase by 1 cap every 3 days to effective dose or max three three times daily",
            "start with 1 capsule at bedtime then increase by 1 capsule a day to effect"
            "inject 10 units subcutaneously every morning and increase by 2 units every other day until fbs 80-120"
            "take one tablet every 12 hours -increase to 2 tablets every 12 hours after 1 week "
            "take 1 tablet by mouth each night at bedtime may increase to 2 tablet"
            "take 6 tablets by mouth on day 1 then decrease by 1 tablet every 3 days"

        One can increment up or down the dose. One can also just set the upper limit: "increase to".

    taper_type              String. Mandatory. What type of taper: 'incremental' 'target'.

    direction               String. Mandatory: 'increase' 'decrease'.

    dose                    DOSE struc. Mandatory. If incremental, then the increment dose (e.g. "increase by 2 caps every day")

    dose_freq               NOT IMPLEMENTED YET. Freq struc. Optional, used only for typ=target. Specifies the dose freq at target level.
                            E.g. "take one tablet every morning can increase to 1 tablet twice daily in 3 days". "twice daily" refers to target freq at target.

    dose_periodicity        Periodicity struc. Optional, used only for typ=target. Specifies the periodicity of dose admin at target level.
                            E.g. "take one tablet every 12 hours -increase to 2 tablets every 12 hours after 1 week as directed"

    increment_periodicity   PERIODICITY struc. Optional. Used only if typ = 'incremental'. Provides the increment value. E.g. "increase by 2 tabs every 3 days"

    offset                  NOT IMPLEMENTED YET. TimeInterval struc: Optional. Used only if typ = 'target' to indicate when one can jump to that target dose.
                            E.g. "take one tablet every 12 hours -increase to 2 tablets every 12 hours after 1 week"

    stop_condition          StopCondition struc. Optional. Specifies when to stop tapering. E.g. "to effect"

    then_flag               True/False. Indicates if there is a Then to the left of Taper, because if there isn't (e.g. 'take 2 tabs daily may increase to 3 daily') we
                            should not have "after that:" between the 2 schedules. But if the previous schedule has Duration and there is "Then" (e.g. "take 2 tabs for 1 week then may increase by 2 tabs"
                            we shloud have "after that:" between the schedules (e.g. "take 2 tabs for 1 week. After that: may increase by ...").

    """

    type_permissible_values = ('incremental', 'target')
    direction_permissible_values = ('increase', 'decrease')

    def __init__(self, taper_type, direction, dose):
        Struc.__init__(self, label = 'TAPER')
        self.taper_type = taper_type
        self.direction = direction
        self.dose = dose
        self.dose_freq = None
        self.dose_periodicity = None
        self.increment_periodicity = None
        self.offset = None
        self.stop_condition = None
        self.then_flag = None


    def is_valid_struc(self):
        return (self.taper_type in self.type_permissible_values and self.direction in self.direction_permissible_values)

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        """

        if self.taper_type != dict_struc.taper_type:
            return None

        if self.direction != dict_struc.direction:
            return None

        dose_score = self.dose.match_dictionary(dict_struc.dose)

        if self.dose_freq and dict_struc.dose_freq:
            dose_freq_score = self.dose_freq.match_dictionary(dict_struc.dose_freq)
        elif self.dose_freq:
            dose_freq_score = 2
        elif dict_struc.dose_freq:
            return None
        else:
            dose_freq_score = None

        if self.dose_periodicity and dict_struc.dose_periodicity:
            dose_periodicity_score = self.dose_periodicity.match_dictionary(dict_struc.dose_periodicity)
        elif self.dose_periodicity:
            dose_periodicity_score = 2
        elif dict_struc.dose_periodicity:
            return None
        else:
            dose_periodicity_score = None

        if self.increment_periodicity and dict_struc.increment_periodicity:
            increment_periodicity_score = self.increment_periodicity.match_dictionary(dict_struc.increment_periodicity)
        elif self.increment_periodicity:
            increment_periodicity_score = 2
        elif dict_struc.increment_periodicity:
            return None
        else:
            increment_periodicity_score = None

        if self.offset and dict_struc.offset:
            offset_score = self.offset.match_dictionary(dict_struc.offset)
        elif self.offset:
            offset_score = 2
        elif dict_struc.offset:
            return None
        else:
            offset_score = None

        if self.stop_condition and dict_struc.stop_condition:
            stop_condition_score = self.stop_condition.match_dictionary(dict_struc.stop_condition)
        elif self.stop_condition:
            stop_condition_score = 2
        elif dict_struc.stop_condition:
            return None
        else:
            stop_condition_score = None

        scores = [dose_score, dose_freq_score, dose_periodicity_score, increment_periodicity_score, offset_score, stop_condition_score]
        return Struc.match_dictionary_calc_multiple_scores(scores)

    def get_numerical_map_to_dict(self, dict_struc):
        """ For numerical variables in the dictionary such as <<NUM_1>> or <<DATE_0>> or <<TIME_0>>, returns a list of pairs
            [(var1, string repr of value1), (var2, value2), ...]. E.g., ('<<NUM_1>>', '1.5'), or ('<<DATE_0>>', '2/3/2012')
        """

        if self.dose and dict_struc.dose:
            dose = self.dose.get_numerical_map_to_dict(dict_struc.dose)
        else:
            dose = []

        if self.dose_freq and dict_struc.dose_freq:
            dose_freq = self.dose_freq.get_numerical_map_to_dict(dict_struc.dose_freq)
        else:
            dose_freq = []

        if self.dose_periodicity and dict_struc.dose_periodicity:
            dose_periodicity = self.dose_periodicity.get_numerical_map_to_dict(dict_struc.dose_periodicity)
        else:
            dose_periodicity = []

        if self.increment_periodicity and dict_struc.increment_periodicity:
            increment_periodicity = self.increment_periodicity.get_numerical_map_to_dict(dict_struc.increment_periodicity)
        else:
            increment_periodicity = []

        if self.offset and dict_struc.offset:
            offset = self.offset.get_numerical_map_to_dict(dict_struc.offset)
        else:
            offset = []

        if self.stop_condition and dict_struc.stop_condition:
            stop_condition = self.stop_condition.get_numerical_map_to_dict(dict_struc.stop_condition)
        else:
            stop_condition = []

        return dose + dose_freq + dose_periodicity + increment_periodicity + offset + stop_condition


    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        direction = self.direction
        taper_type = '_by_' if self.taper_type == 'incremental' else '_to_'
        dose = self.dose.string_key()

        return direction + taper_type + dose

    def is_semantically_incompatible_with_given_sem(self, schedule):
        """ Returns True always (if there is a dose in any of the events of the schedule) because we want to start a new Schedule whenever there is a Taper instruction.

        """

        for event in schedule.events:
            if event.dose:
                return True
        else:
            return False

    def quick_print_struc(self):

        direction = self.direction
        taper_type = '_by_' if self.taper_type == 'incremental' else '_to_'
        dose = self.dose.quick_print_struc()
        dose_freq = ('/dose_freq-' + self.dose_freq.quick_print_struc()) if self.dose_freq else ''
        dose_periodicity = ('/dose_period-' + self.dose_periodicity.quick_print_struc()) if self.dose_periodicity else ''
        increment_periodicity = ('/increment_interval-' + self.increment_periodicity.quick_print_struc()) if self.increment_periodicity else ''
        offset = ('/after-' + self.offset.quick_print_struc()) if self.offset else ''
        stop_condition = ('/stop-' + self.stop_condition.quick_print_struc()) if self.stop_condition else ''

        return direction + taper_type + dose + dose_freq + dose_periodicity + increment_periodicity + offset + stop_condition



class Route(Struc):
    """ Represents route, e.g by mouth, intravenously, etc.
    value                       String for route, capitalized

    permissible_values          List
    value_to_likely_directive   Dictionary mapping route.value to possible directive values (as txt).
    value_to_likely_form        Dictionary mapping route.value to possible Form values (as txt).
    """

    permissible_values = ('externally', 'intramuscularly', 'intranasally', 'orally', 'rectally', 'subcutaneously', 'sublingually',
                          'topically', 'transdermally', 'vaginally')

    value_to_likely_directive = {'externally': 'apply',
                                 'intramuscularly': 'inject',
                                 'intranasally': 'use',
                                 'orally': 'take',
                                 'rectally': 'apply',
                                 'subcutaneously': 'inject',
                                 'sublingually': 'dissolve',
                                 'topically': 'apply',
                                 'vaginally': 'insert'}

    value_to_likely_form = {'externally': 'patch',
                            'intranasally': 'spray',
                            'orally': 'tablet',
                            'rectally': 'suppository',
                            'subcutaneously': 'unit',
                            'sublingually': 'tablet',
                            'topically': 'patch',
                            'vaginally': 'applicatorful'}

    def __init__(self, value, constituents):
        Struc.__init__(self, label = 'ROUTE', constituents = constituents)
        self.value = value

    def is_valid_struc(self):
        test_result = self.value in self.permissible_values
        return test_result

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        """

        if self.value == dict_struc.value:
            return 4
        else:
            return 0

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """


        if self.value == 'orally':
            # because with tablets "by mouth" is redundant enough, its absence in the sig should not disquality an atom with "by mouth" from being matched to the sig.
            value = ''
        else:
            value = self.value

        return value

    def quick_print_struc(self):
        return self.value

class Site(Struc):
    """ Represents the site, e.g "apply to affected areas", "use 2 sprays to each nostril"
    value       string for site
    relation    string (optional), a preposition of how to do the thing to the site, (to, into, on, in)
    which       string (optional). Used when the site consists of a potential multiplicity of specific sites,
                e.g. left nostril, each nostril, into one nostril, affected eye, each ear.

    permissible_values  list for permissible values for site name
    which_permissible_values    left|right|affected|both     (as in "affected eye", not "affected area" which is a site value)
    value_to_likely_directive   Dict that maps site.value to typical directives using them. It is used to deduce the directive when the site is given and directive is not or is non-standard (not found in dictionary)
    """
    permissible_values = ('affected area', 'acne areas', 'back', 'body', 'cheek', 'diaper area', 'ear', 'eye', 'eyelid', 'eyelash', 'face', 'feet', 'foot',
                          'hair', 'hands', 'joint', 'knees', 'legs', 'lesion', 'lower back', 'mouth', 'nails', 'nostril', 'painful area', 'rash', 'scalp',
                          'skin', 'vagina', 'warts', 'wound')

    which_permissible_values = ('left', 'right', 'each', 'affected', 'upper', 'lower', 'one')       # upper/lower pertain to eyelids
    value_to_likely_directive = {'acne areas': 'apply',
                                 'affected area': 'apply',
                                 'back': 'apply',
                                 'body': 'apply',
                                 'cheek': 'apply',
                                 'diaper area': 'apply',
                                 'ear': 'instill',
                                 'eye': 'instill',
                                 'eyelid': 'apply',
                                 'eyelash': 'apply',
                                 'face': 'apply',
                                 'feet': 'apply',
                                 'foot': 'apply',
                                 'hair': 'apply',
                                 'hands': 'apply',
                                 'joint': 'apply',
                                 'knees': 'apply',
                                 'legs': 'apply',
                                 'lesion': 'apply',
                                 'lower back': 'apply',
                                 'nails': 'apply',
                                 'nostril': 'spray',
                                 'painful area': 'apply',
                                 'rash': 'apply',
                                 'scalp': 'apply',
                                 'skin': 'apply',
                                 'vagina': 'apply',
                                 'warts': 'apply',
                                 'wound': 'apply'}

    def __init__(self, value, relation, constituents, which = None):
        Struc.__init__(self, label = 'SITE', constituents = constituents)
        self.value = value
        self.relation = relation
        self.which = which

    def is_valid_struc(self):
        test_result = self.value in self.permissible_values and (self.which is None or self.which in self.which_permissible_values )
        return test_result

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        0   unacceptable
        1   not right, but acceptable if no other choices
        2   acceptable but defective
        3   acceptable but not perfect
        4   perfect match
        """

        if self.value == dict_struc.value:
            value_score = 4
        elif dict_struc.value == 'affected area':
            value_score = 1.5
        else:
            return 0

        if self.which == dict_struc.which:
            which_score = 4
        elif self.which in ('left', 'right') and dict_struc.which == 'affected':
            which_score = 1
        else:
            return 0


        scores = [value_score, which_score]
        return Struc.match_dictionary_calc_multiple_scores(scores)



    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        if self.value == 'affected area':
            value = ''
        else:
            value = self.value

        return value

    def quick_print_struc(self):
        which = (self.which + '_' ) if self.which else ''
        return self.relation + '_' + which +  self.value

class Timing(Struc):
    """ Represents the timing of action, i.e. when it occurs: e.g. before meals, 2 hours after breakfast, 1 hour before sex, etc.

    typ         String. Type of Timing (actually type of landmark).
                    - time_of_day   landmark is of the form 'HH:MM' in military time
                    - food_meals    main_meals | food_landmarks - a main meal ('breakfast', 'lunch', 'dinner') or food ('food', 'meal', 'food or milk', 'empty stomach')
                    - day_part      dayparts | bedtime:     'morning', 'afternoon', 'evening', 'night', 'bedtime'
                    - event         med_events | misc_landmarks: 'sex', 'dental appointment', 'appointment', 'diaper change', 'bowel_movement', 'loose_stool'

    landmark    string representing the event with respect to which the action
                needs to take place, e.g. night, breakfast, meal, empty stomach, food and drink, appointment

                Optional:
    relation    (optional) String. Typically a proposition, repr. relationship to the landmark, e.g. with, after, during, on, before
                Used only for relationships with discrete events (e.g. meals), not day parts (e.g. "evening").

    offset      (optional) TIMEINTERVAL Struc representing the number of time units offset from landmark,
                eg. "2 hours before breakfast"

    start_flag  (optional) None/True flag representing whether the Timing is the start of the series of administration or not.
                E.g. "instill 2 drops in the right eye starting 2 days before surgery."

    every_flag  (optional) None/True flag representing that the action should take place at the specified
                landmark at each of it's occurences. E.g. "with each meal and at bedtime", "every night at bedtime"
                If the landmark is a once-a-day event (breakfast, bedtime, etc) then the presence of "every"
                and "each" implies that the action should happen every day, i.e. we need to add
                something like "Take this medicine every day" if it is not expressed elsewhere in the Sig.
                But if the landmark can happen several times daily, e.g. "after every meal", then the daily freq
                is possibly different from 1.

    landmark_permissible_values     Set of strings, used for validation and in self.incompatible_landmarks() function
    relation_permissible_values     Set of strings
    """

    main_meals = set(['breakfast', 'lunch', 'dinner'])
    food_landmarks = set(['empty stomach', 'food or milk', 'food', 'meal'])
    dayparts = set(['morning', 'afternoon', 'evening', 'night'])
    sleep_events = set(['bedtime', 'awakening'])
    med_events = set(['dental appointment', 'appointment', 'procedure'])
    misc_landmarks = set(['sex', 'onset_of_headache', 'diaper change', 'bowel_movement', 'loose_stool'])

    landmark_permissible_values = main_meals | food_landmarks | dayparts | sleep_events | med_events | misc_landmarks
    relation_permissible_values = ('at', 'before', 'after')

    pattern_time_of_day = re.compile('\d\d:\d\d|<<TIME\_\d+>>')   # military time 'HH:MM' or Quant of num_type = "var" and var_type = 'time' (e .g <<TIME_0>>)


    def __init__(self, landmark, relation, typ = None, offset = None, start_flag = None, every_flag = None):
        Struc.__init__(self, label = 'TIMING')
        self.typ = typ
        self.landmark = landmark
        self.relation = relation
        self.offset = offset
        self.start_flag = start_flag
        self.every_flag = every_flag

        if not typ:
            if landmark in Timing.main_meals | Timing.food_landmarks:
                self.typ = 'food_meals'
            elif landmark in Timing.dayparts | Timing.sleep_events:
                self.typ = 'day_part'
            elif landmark in Timing.med_events | Timing.misc_landmarks:
                self.typ = 'event'
            elif Timing.pattern_time_of_day.match(landmark):
                self.typ = 'time_of_day'

        constituents = [landmark]

        if relation:
            constituents += [relation]
        if offset:
            constituents += [offset]
        if every_flag:
            constituents += ['every']

        self.constituents = constituents

        if debug and not self.is_valid_struc():
            raise Exception('Landmark %s is not on the list of landmark_permissible_values' % landmark)

    def is_valid_struc(self):
        landmark_ok = self.landmark in self.landmark_permissible_values or Timing.pattern_time_of_day.match(self.landmark)
        relation_ok = not self.relation or self.relation in self.relation_permissible_values
        return landmark_ok and relation_ok



    @staticmethod
    def incompatible_landmarks(landmark1, landmark2):
        """ Returns T/F. Two landmarks are incompatible if they can't be part of the same AdminEvent.

        E.g., morning and bedtime. But night and bedtime are compatible. Even 'breakfast" and 'empty stomach' are compatible
        because one can say "on an empty stomach before breakfast".  """

        if landmark1 == landmark2:
            return False

        landmarks = set((landmark1, landmark2))


        if landmarks <= (Timing.main_meals | Timing.sleep_events):
            return True
        elif landmarks <= Timing.dayparts:
            return True
        elif 'bedtime' in landmarks and set(['morning', 'afternoon']) & landmarks:
            return True
        elif 'awakening' in landmarks and set(['afternoon', 'evening', 'night']) & landmarks:
            return True
        elif 'breakfast' in landmarks and set(['afternoon', 'evening', 'night']) & landmarks:
            return True
        elif 'lunch' in landmarks and set(['morning', 'evening', 'night']) & landmarks:
            return True
        elif 'dinner' in landmarks and set(['morning', 'afternoon']) & landmarks:
            return True
        elif 'empty stomach' in landmarks and set(['food or milk', 'food', 'meal']) & landmarks:
            return True
        elif Timing.pattern_time_of_day.match(landmark1) and Timing.pattern_time_of_day.match(landmark2):
            return True
        else:
            return False

    def is_semantically_incompatible_with_given_sem(self, admin_event):
        """ Returns True if the Timing structure is incompatible with at least one of the Timing structures on the admin_event.timing (which is a list
            of Timing strucs). I.e. Self is incompatible with the timing list of an AdminEvent if Self can't be added to that AdminEvent.timing -- it signals
            a new AdminEvent.

        Two Timing structures are incompatible if they both have landmarks and those landmarks are incompatible (e.g. "morning" vs. "evening")
        or they both have landmarks and these landmarks are identical, but have different offsets.
        If either of them is missing a landmark, they are compatible.
        """

        landmark = self.landmark
        offset = self.offset

        if not landmark:
            return False

        timing_list = admin_event.timing

        for timing in timing_list:
            current_landmark = timing.landmark
            current_offset = timing.offset
            if not current_landmark:
                continue
            if self.incompatible_landmarks(current_landmark, landmark):
                return True
            if current_landmark == landmark and current_offset and current_offset.quant and offset and offset.quant and current_offset.quant != offset.quant:
                return True
        return False

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        0   unacceptable
        1   not right, but acceptable if no other choices
        2   acceptable but defective
        3   acceptable but not perfect
        4   perfect match
        """


        if self.landmark == dict_struc.landmark:
            landmark_score = 4
        elif self.typ == 'time_of_day' and dict_struc.typ == 'time_of_day' and 'TIME' in dict_struc.landmark:
            landmark_score = 4
        elif self.landmark == 'dental appointment' and dict_struc.landmark == 'appointment':
            landmark_score = 2
        else:
            return 0

        if self.relation == dict_struc.relation:
            relation_score = 4
        elif self.landmark in Timing.sleep_events | Timing.misc_landmarks:
            relation_score = 3
        else:
            return 0

        if self.offset and not dict_struc.offset:
            return 0
        if dict_struc.offset and not self.offset:
            return 0
        if self.offset and dict_struc.offset:
            offset_score = self.offset.match_dictionary(dict_struc.offset)
        else:
            offset_score = 4

        if self.every_flag == dict_struc.every_flag:
            every_flag_score = 4
        elif not self.every_flag and not dict_struc.every_flag:
            # sometimes one flag is None and the other is False
            every_flag_score = 4
        else:
            every_flag_score = 2

        if self.start_flag == dict_struc.start_flag:
            start_flag_score = 4
        elif not self.start_flag and not dict_struc.start_flag:
            # sometimes one flag is None and the other is False
            start_flag_score = 4
        else:
            start_flag_score = 1


        scores = [landmark_score, relation_score, offset_score, start_flag_score, every_flag_score]
        return Struc.match_dictionary_calc_multiple_scores(scores)

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        if self.offset:
            offset = '-' + self.offset.string_key()
        else:
            offset = ''

        return self.typ + offset


    def get_numerical_map_to_dict(self, dict_struc):
        """ For numerical variables in the dictionary such as <<NUM_1>> or <<DATE_0>> or <<TIME_0>>, returns a list of pairs
            [(var1, string repr of value1), (var2, value2), ...]. E.g., ('<<NUM_1>>', '1.5'), or ('<<DATE_0>>', '2/3/2012')
        """

        if self.offset:
            offset = self.offset.get_numerical_map_to_dict(dict_struc.offset)
        else:
            offset = []

        if self.typ == 'time_of_day' and dict_struc.typ == 'time_of_day' and 'TIME' in dict_struc.landmark:
            variable = dict_struc.landmark # landmark looks like <<TIME_N>>
            value = self.landmark
            landmark = [(variable, value)]
        else:
            landmark = []

        return landmark + offset

    def quick_print_struc(self):
        offset = (self.offset.quick_print_struc() + '-') if self.offset else ''
        relation = (self.relation + '_') if self.relation else ''
        typ = 'Type ' + self.typ + ':'
        every = 'every_' if self.every_flag else ''
        start = 'starting_' if self.start_flag else ''

        return typ + start + offset + relation + every + self.landmark

class Calendar_Event(Struc):
    """ Represents a specific, calendar-tied event, e.g. now, tomorrow, on day 2, on days 2-7, etc. See Specific_Day for days of the week or specific dates.

    typ     String: types are:
                    'numeric'   e.g. "first day" or "days 2 - 5"
                    'relative'  e.g. 'now', 'today', 'tomorrow'


            Optional:
    quant       (optional) QUANT. Occurs only when typ == 'numeric'
    time_unit   (optional) TIMEUNIT. Occurs only when typ == 'numeric'. E.g. = 'day' if event is 'days 2-7'
    value       (optional) String. Has value if typ != 'numeric'
                If typ == 'relative': value is in ('now', 'tomorrow', 'today').

    """

    typ_permissible_values = ('numeric', 'relative')

    def __init__(self, typ):
        Struc.__init__(self, label = 'CALENDAR_EVENT')
        self.typ = typ
        self.constituents = []
        self.quant = None
        self.time_unit = None
        self.value = None


    def is_valid_struc(self):
        return self.typ in self.typ_permissible_values


    @staticmethod
    def incompatible_calendar_event(calendar_event1, calendar_event2):
        """ Returns T/F. Two calendar_events are incompatible if they can't be part of the same Schedule.

        E.g., "now" and "tomorrow" are incompatible.
        Two calendar_event structures are incompatible if they are of the same type but not identical or if one is numeric and another is relative.
        """

        if calendar_event1.typ == calendar_event2.typ:
            if calendar_event1.pointwise_equal(calendar_event2):
                return False
            else:
                return True

        typs = set((calendar_event1.typ, calendar_event2.typ))
        if typs == set(['numeric', 'relative']):
            return True
        else:
            return False

    def is_semantically_incompatible_with_given_sem(self, schedule):
        """ Returns True if the Calendar_Event structure is incompatible with at least one of the Calendar_Event structures on the schedule.calendar_event
            (which is a list of Calendar_Event strucs). I.e. Self is incompatible with the calendar_event list of Schedule if Self can't be added to that
            Schedule.calendar_event -- it signals a new Schedule.

        Two calendar_event structures are incompatible if they are of the same type but not identical or if one is numeric and another is relative.
        """

        calendar_event_list = schedule.calendar_event

        for calendar_event in calendar_event_list:
            current_typ = calendar_event.typ
            if not current_typ:
                continue
            if self.incompatible_calendar_event(self, calendar_event):
                return True
        return False

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        0   unacceptable
        1   not right, but acceptable if no other choices
        2   acceptable but defective
        3   acceptable but not perfect
        4   perfect match
        """

        if self.typ == dict_struc.typ:
            typ_score = 4
        else:
            return 0

        if self.quant:
            quant_score = self.quant.match_dictionary(dict_struc.quant)
        else:
            quant_score = 4

        if self.time_unit:
            time_unit_score = self.time_unit.match_dictionary(dict_struc.time_unit)
        else:
            time_unit_score = 4

        if self.value == dict_struc.value:
            value_unit_score = 4
        else:
            return 0

        scores = [typ_score, quant_score, time_unit_score, value_unit_score]
        return Struc.match_dictionary_calc_multiple_scores(scores)


    def get_numerical_map_to_dict(self, dict_struc):
        """ For numerical variables in the dictionary such as <<NUM_1>> or <<DATE_0>> or <<TIME_0>>, returns a list of pairs
            [(var1, string repr of value1), (var2, value2), ...]. E.g., ('<<NUM_1>>', '1.5'), or ('<<DATE_0>>', '2/3/2012')
        """

        if self.quant:
            quant = self.quant.get_numerical_map_to_dict(dict_struc.quant)
            return quant
        else:
            return []

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        if self.quant:
            quant = '-' + self.quant.string_key()
        else:
            quant = ''

        if self.time_unit:
            time_unit = '-' + self.time_unit.string_key()
        else:
            time_unit = ''

        return self.typ + quant + time_unit


    def quick_print_struc(self):
        numeric = (self.time_unit.quick_print_struc() + '-' + self.quant.quick_print_struc()) if self.typ == 'numeric' else ''
        value = self.value if self.typ != 'numeric' else ''
        return numeric + value


class Specific_Day(Struc):
    """ Represents a specific day or day range: Saturday, Weekend, Mon-Friday, 2/3/2013

    typ     String: types are:
                    'day_of_week' i.e. Monday - Sunday or range: weenend, weekdays
                    'date'  Date strtucture

            Optional:
    day_of_week     (optional) String. If typ == 'day_of_week': value = name of day of the week ('monday',... 'sunday') or "weekend" or "weekdays".
    date            (optional) DATE struc if typ == 'date'

    """

    typ_permissible_values = ('day_of_week', 'date')
    day_of_week_permissible_values = ('monday', 'tuesday', 'wednesday', 'thursday', 'friday', 'saturday', 'sunday', 'weekdays', 'weekend')

    def __init__(self, typ, day_of_week = None):
        Struc.__init__(self, label = 'SPECIFIC_DAY')
        self.typ = typ
        self.day_of_week = day_of_week
        self.date = None
        self.constituents = []

    def is_valid_struc(self):
        if self.typ not in self.typ_permissible_values:
            return False
        if self.day_of_week and self.day_of_week not in self.day_of_week_permissible_values:
            return False
        return True


    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        0   unacceptable
        1   not right, but acceptable if no other choices
        2   acceptable but defective
        3   acceptable but not perfect
        4   perfect match
        """

        if self.typ == dict_struc.typ:
            typ_score = 4
        else:
            return 0

        if self.day_of_week == dict_struc.day_of_week:
            day_of_week_score = 4
        else:
            return None

        if self.date and dict_struc.date:
            date_score = self.date.match_dictionary(dict_struc.date)
        elif self.date or dict_struc.date:
            return None
        else:
            date_score = 4


        scores = [typ_score, day_of_week_score, date_score]
        return Struc.match_dictionary_calc_multiple_scores(scores)


    def get_numerical_map_to_dict(self, dict_struc):
        """ For numerical variables in the dictionary such as <<NUM_1>> or <<DATE_0>> or <<TIME_0>>, returns a list of pairs
            [(var1, string repr of value1), (var2, value2), ...]. E.g., ('<<NUM_1>>', '1.5'), or ('<<DATE_0>>', '2/3/2012')
        """

        if self.date:
            date = self.date.get_numerical_map_to_dict(dict_struc.date)
            return date
        else:
            return []

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        if self.day_of_week:
            day_of_week = '-' + self.day_of_week
        else:
            day_of_week = ''

        return self.typ + day_of_week


    def quick_print_struc(self):
        day_of_week = self.day_of_week if self.day_of_week else ''
        date = self.date.quick_print_struc() if self.date else ''
        return day_of_week + date



class Indication(Struc):
    """ For what condition is the drug prescribed, e.g. "severe abdominal pain", "to reduce pain or fever"

    condition       String  e.g. "pain", "allergy"
    """

    permissible_values = ('acid reflux', 'acne', 'agitation', 'allergies', 'anemia', 'anxiety', 'arthritis', 'asthma', 'bladder', 'blood circulation', 'blood thinner', 'blood pressure',
                          'bones', 'breathing', 'calcium', 'cholesterol', 'concentration', 'congestion',
                          'constipation', 'cough', 'cough and congestion', 'cough/wheeze', 'cramps', 'dementia', 'depression',
                          'diabetes', 'dialysis', 'diarrhea', 'digestion', 'dizziness', 'dry eyes', 'dry skin',
                          'dryness', 'eczema', 'edema', 'fever', 'fungal infection', 'gas', 'glaucoma', 'gout', 'headache', 'heartburn', 'heart rate',
                          'heart', 'hemorrhoids', 'hiccups', 'high blood pressure', 'incontinence', 'indigestion', 'infection', 'inflammation', 'iron', 'irritation', 'itching', 'leg swelling',
                          'leg', 'lungs', 'memory', 'migraine', 'mood', 'muscle spasms', 'nausea', 'nausea and vomiting',
                          'nervousness', 'osteoporosis', 'pain', 'pain or fever', 'panic attack', 'potassium', 'prostate', 'rash', 'rhinitis', 'runny nose',
                          'seizures', 'shortness of breath', 'sleep/insomnia', 'sore throat', 'spasms', 'stomach', 'stomach cramps', 'stuffy nose',
                          'swelling', 'thyroid', 'tremors', 'triglycerides', 'urinary tract infection', 'vertigo', 'vomiting', 'wheesing cough', 'wheezing')

    def __init__(self, condition, constituents):
        Struc.__init__(self, label = 'INDICATION')
        self.condition = condition
        self.constituents = [constituents]

    def is_valid_struc(self):
        test_result = self.condition in self.permissible_values
        return test_result

    def is_semantically_incompatible_with_given_sem(self, admin_event):
        """ Returns True if the Indication structure is incompatible with at least one of the Indication structures on the Instruction.indication list
            (which is a list of Indication strucs). I.e. Self is incompatible with the indication list of an Instruction if Self can't be added to that Instruction.indication -- it signals
            a new Instruction.

        For now, this is a dummy returning False.
        """

        return False


    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        """

        if self.condition == dict_struc.condition:
            return 4
        else:
            return 0

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        return self.condition

    def quick_print_struc(self):
        return 'for-' + self.condition


class IndicationPain(Indication):
    """ For what type of pain is this prescribed.

    pain_type       String: typically location of pain, e.g. "stomach", "chest", etc.
    severity        True/False (e.g "severe stomach pain")
    verb            String: e.g. "reduce"
    permissible_values  List

    """

    permissible_values = ('back', 'breakthrough', 'chest', 'ear', 'joint', 'knee', 'nerve', 'stomach', '')

    def __init__(self, pain_type, severity, verb, constituents):
        Indication.__init__(self, condition = 'pain', constituents = constituents)
        self.pain_type = pain_type
        self.severity = severity
        self.verb = verb

    def is_valid_struc(self):
        test_result = self.pain_type in self.permissible_values
        return test_result


    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        0   unacceptable
        1   not right, but acceptable if no other choices
        2   acceptable but defective
        3   acceptable but not perfect
        4   perfect match
        """

        if dict_struc.condition != 'pain':
            return 0
        elif self.pain_type == dict_struc.pain_type:
            pain_type_score = 4
        elif not dict_struc.pain_type:
            pain_type_score = 1
        else:
            return 0

        if self.severity == dict_struc.severity:
            severity_score = 4
        elif not dict_struc.severity:
            severity_score = 2
        else:
            severity_score = 1

        if not self.severity and not dict_struc.severity:
            # We don't want the absence of "severety" to tamper down the difference between a really good match of pain types between different dictionary entries.
            return pain_type_score
        else:
            scores = [pain_type_score, severity_score]
            return Struc.match_dictionary_calc_multiple_scores(scores)


    def quick_print_struc(self):
        verb = (self.verb + '-') if self.verb else ''
        severity = 'severe-' if self.severity else ''
        pain_type = (self.pain_type + '-') if self.pain_type else ''
        return verb + severity + pain_type + 'pain'


class AndConj(Struc):
    """ Represents conjunction "and"

    E.g. conjunction of TIMING: "take 2 pills in the morning and 2 in the evening"
    conjunction of Indicantions: "take for heart and for blood pressure."

    """

    def __init__(self, constituents):
        Struc.__init__(self, label = 'AND_CONJ', constituents = constituents)

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        """
        return 4

    def quick_print_struc(self):
        return 'and'


class ThenChrono(Struc):
    """ Represents chronological connective "then" (e.g. "take 2 tabls days 1-5, then 1 thereafter").
    """

    def __init__(self, constituents):
        Struc.__init__(self, label = 'THEN_CHRONO', constituents = constituents)

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        """
        return 4

    def quick_print_struc(self):
        return 'then'


class Anaphora(Struc):
    """ Represents anaphoric reference to a prior concept. E.g. "Take this medicine in the evenings" or "this is your first dose" or "do this <<NUM_0>> hours prior to appointment"
        Used almost exclusively to process Dictionary.

    value           this_medicine | this_dose | do_this | this_is | this_will

    """

    permissible_values = ('this_medicine', 'this_dose', 'do_this', 'this_is', 'this_will')

    def __init__(self, value, constituents):
        Struc.__init__(self, label = 'ANAPHORA', constituents = constituents)
        self.value = value


    def is_valid_struc(self):
        return self.value in self.permissible_values

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        """

        if self.value == dict_struc.value:
            return 4
        else:
            return 3.9  # Because we tend to insert into the Event of the sig only 1 Anaphora ("This medicine") but in reality the use of "do_this" and "this_is" will take care of its own matching.

    def quick_print_struc(self):
        return 'Anaphora:' + self.value


class AsDirected(Struc):
    """ Represents directive modifiers such as "as directed", "per instructions of md", "follow packaging instructions"

    verb            string. The action of the prescribing authority: E.g. directed, instructed, explained, prescribed
    authority       string. E.g. doctor, package, etc.
    exactly_flag    True/False if "exactly" is present, e.g. "use exactly as prescribed"
    your_flag       True/False if "your doctor" is present.
    or_flag         True/False if "or as directed" type of struc is present.
    directive       (Optional) String: The action you are to perform with the prescription (similar to DIRECTIVE): take, use, give, inject, etc.

    The difference between verb, directive, and authority is illustrated here: "take as directed by your doctor":
    verb = 'directed'
    authority = 'doctor'
    directive = 'take'

    """

    def __init__(self, verb, constituents, authority = None, exactly_flag = None, your_flag = None, or_flag = None):
        Struc.__init__(self, label = 'AS_DIRECTED', constituents = constituents)
        self.verb = verb
        self.authority = authority
        self.exactly_flag = exactly_flag
        self.your_flag = your_flag
        self.or_flag = or_flag
        self.directive = None

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        0   unacceptable
        1   not right, but acceptable if no other choices
        2   acceptable but defective
        3   acceptable but not perfect
        4   perfect match
        """

        if self.verb == dict_struc.verb:
            verb_score = 4
        else:
            verb_score = 2

        if self.authority == dict_struc.authority:
            authority_score = 4
        elif not dict_struc.authority:
            authority_score = 2
        elif not self.authority and dict_struc.authority == 'doctor':   # if sig says "take as directed" and dict says "take as directed by doctor", it is a bearable match,
                                                                        # especially if there is a match on a rare aspect ssuch as "or_flag"
            authority_score = 1.5
        else:
            return 0

        if self.exactly_flag == dict_struc.exactly_flag:
            exactly_score = 4
        else:
            exactly_score = 2

        if self.your_flag == dict_struc.your_flag:
            your_score = 4
        else:
            your_score = 3

        if self.or_flag == dict_struc.or_flag:
            or_score = 4
        elif not dict_struc.or_flag:
            or_score = 2
        else:
            return 0

        if self.directive == dict_struc.directive:
            directive_score = 4
        elif not dict_struc.directive or dict_struc.directive == 'use':
            directive_score = 3
        else:
            return 0

        scores = [verb_score, authority_score, exactly_score, your_score, or_score, directive_score]
        return Struc.match_dictionary_calc_multiple_scores(scores)


    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        authority = self.authority
        if authority == 'doctor':
            authority = ''

        or_flag = 'OR' if self.or_flag else ''
        value = '-'.join([prop for prop in (or_flag, authority) if prop])

        value = value if value else True
        return value

    def quick_print_struc(self):
        authority = ('_by_' + ('your_' if self.your_flag else '') + self.authority ) if self.authority else ''
        exactly = 'exactly-' if self.exactly_flag else ''
        or_disjunction = 'or-' if self.or_flag else ''
        directive = (self.directive + '_') if self.directive else ''
        verb = ('as_' + self.verb) if self.verb else ''
        return or_disjunction + directive + exactly + verb + authority

class AsNeeded(Struc):
    """ Represents directive modifiers such as "as needed", "only if needed"

    only_flag   True/False if "only" is present, e.g. "use only when needed"
    if_flag     True/False if "if" is present, e.g. "use if you need it"
    """

    def __init__(self, constituents, only_flag = None):
        Struc.__init__(self, label = 'AS_NEEDED', constituents = constituents)
        self.only_flag = only_flag
        self.if_flag = False


    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        """

        if self.only_flag == dict_struc.only_flag:
            only_score = 4
        else:
            only_score = 1

        if self.if_flag == dict_struc.if_flag:
            if_score = 4
        else:
            if_score = 2

        scores = [only_score, if_score]
        return Struc.match_dictionary_calc_multiple_scores(scores)

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        return True

    def quick_print_struc(self):
        if_flag = 'if_' if self.if_flag else ''
        only_flag = 'only_' if self.only_flag else ''
        return if_flag + only_flag + 'as_needed'


class SeeInstructions(Struc):
    """ Represents instructions such as "see attached label for directions" or "see note"

    authority       string. What the person is to see, e.g. label, directions, note. Can be empty (e.g. "See attached."
    attached_flag   True/False if "attached" is present, e.g. "see attached note"

    """

    def __init__(self, authority, constituents, attached_flag = None):
        Struc.__init__(self, label = 'SEE_INSTRUCTIONS', constituents = constituents)
        self.authority = authority
        self.attached_flag = attached_flag

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        """

        if self.authority == dict_struc.authority:
            authority_score = 4
        elif not dict_struc.authority:
            authority_score = 2
        else:
            return 0

        if self.attached_flag == dict_struc.attached_flag:
            attached_score = 4
        elif not dict_struc.attached_flag:
            attached_score = 3
        else:
            attached_score = 1

        scores = [authority_score, attached_score]
        return Struc.match_dictionary_calc_multiple_scores(scores)


    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        return True

    def quick_print_struc(self):
        attached = 'attached_' if self.attached_flag else ''
        authority = self.authority
        return 'see_' + attached + authority



class StopCondition(Struc):
    """ The condition which has to trigger discontinuation of using the medication or stopping tapering (dose increase/decrease)

    condition:          String: e.g "until gone", or "to effect"

    """

    permissible_values = ('if_lightheaded', 'if_side_affects_appear', 'to_effect', 'unless_new_symptoms_appear', 'unless_symptoms_worsen',
                          'until_gone', 'until_return_from_trip', 'until_symptoms_relieved', 'until_therapy_completed')

    def __init__(self, condition, constituents):
        Struc.__init__(self, label = 'STOP_CONDITION', constituents = constituents)
        self.condition = condition

    def is_valid_struc(self):
        return self.condition in self.permissible_values

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        """

        if self.condition == dict_struc.condition:
            return 4
        else:
            return 0

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        return self.condition

    def quick_print_struc(self):
        return self.condition

class DiscardRemanider(Struc):
    """ Represents instructions such as "discrd the rest"
    """

    def __init__(self, constituents):
        Struc.__init__(self, label = 'DISCARD_REMAINDER', constituents = constituents)

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        """
        return 4

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        return True

    def quick_print_struc(self):
        return 'discard_remainder'

class Miscellaneous(Struc):
    """ Represents instructions such of several types.

    typ                     String in permissible_value list
    permissible_values      tuple of permissible values for typ

            Examples: for typ = '911': e.g. "If no relief call 911"
    """

    permissible_values = ('call_911', 'call_911_if_condition_persists', 'administer_in_office', 'administer_in_pharmacy', 'inject_in_office', 'inject_at_pharmacy', 'for_insulin_injections',
                          'antibiotic', 'do_not_swallow', 'rinse_mouth_after_use', 'vitamin', 'diuretic', 'supplement')

    def __init__(self, typ, constituents):
        Struc.__init__(self, label = 'MISCELLANEOUS', constituents = constituents)
        self.typ = typ

    def is_valid_struc(self):
        test_result = self.typ in self.permissible_values
        return test_result

    def match_dictionary(self, dict_struc):
        """ Returns a number 0-4 indicating the quality of the match.
        """

        if self.typ == dict_struc.typ:
            return 4
        else:
            return 0

    def string_key(self):
        """ Returns a string that represents the essential elements of the struc that are minimally required to be matched to a dictionary atom.

        Used in get_key_set(struc_list). The idea is that all atoms compatible with the sig should have a set of key_set that is a subset of the keys
        of the Sig.

        If returns True, then the actual key would be the label.
        If returns None or '', there is no key
        For all other returned values, get_key() will use struc.label + '/' + key
        """

        return self.typ

    def quick_print_struc(self):
        return self.typ


class RuleI(object):
    """ Rule Interface object. Needs to be subclassed.

    Rules for tranforming sig's parse structures to extract information and translform parses.

    nume                            STR: Name of individual rule instance, e.g. "timing".

    rule_type                       STR: Rule type, typically the name of the rule class (e.g. 'Rule_ExtractStrucs')

    rule_properties                 SET: set of rule properties (as strings), to enable to extract various subsets of rules as  needed.
                                    e.g. rule_type includes possible properties such as 'preprocessing', 'struc_identification'.

    change_parse_in_place           TRUE/FALSE: If True, the procedure does not create new parses, but modifies in place of existing parses.
                                    If False, then if the rule applies, a new parse is created with the results of that rule application, while the old rule is retained.

"""

    def __init__(self, name, change_parse_in_place):
        self.name = name
        self.rule_type = ''
        self.rule_properties = set()
        self.change_parse_in_place = change_parse_in_place

    def apply_rule(self, sig):
        pass

class Rule_ExtractStrucs(RuleI):
    """ These rules extract structure from raw text.

    search_patterns                 A list of reg expr search patterns to be used in the search_proc.
                                    We sometimes use several related search patterns to be used with one replacement_proc if the
                                    patterns differ but they look for the same thing and use similar groupdict.

    search_proc                     A search procedure (typically re. search) that are applied to flattened parse, or custom procedures that
                                    behave similarly (e.g. extraction and evaluation of numerics via external module).

    replacement_proc                A procedure that takes as parameters:
                                    a. search pattern's search_object
                                    b. left context string
                                    c. right context string
                                    d. parse (this is sometimes needed for scanning the rest of the parse strucs while validating the rule application)
                                    and either returns None or replaces some subset of Sig parses' strucs with new strucs
    """

    def __init__(self, name, search_patterns, search_proc, replacement_proc, rule_properties, change_parse_in_place = None):

        RuleI.__init__(self, name, change_parse_in_place)
        self.search_patterns = search_patterns
        self.search_proc = search_proc
        self.replacement_proc = replacement_proc
        self.rule_type = 'Rule_ExtractStrucs'
        self.rule_properties |= rule_properties

    def apply_rule(self, sig):
        """
        Applies the rule by trying to match the search patterns to each parse struc one by one, and for each match
        transforms that struc by applying self.replacement_proc to the old struc.

        self.replacement_proc takes 3 arguments: the found_obj, left and right contexts. If it returns None the replacement is abandoned.

        If the flag self.change_parse_in_place is True, the affected parse is changed in place (e.g. for the cases when we are sure about the
        unambiguous nature of the match and we don't want to create too many alternative parses). Typically this happens in pre-processing.

        If the flag self.change_parse_in_place is False, each time there is a match of a struc, we leave the old parse alone and create a
        new parse with the old struc replaced by new transformation. The new parse(s) are added to the parse queue.

        Whenever a parse is created or changed, we check to see if it duplicates an existing parse. This is important to make sure that
        order of rule application for independent rules doesn't create proliferation of redundant parses.

        """

        new_parses = []
        for parse_num, parse in enumerate(sig.parses):
            parse_changed_or_new_parse_found = False
            for pattern_num, pattern in enumerate(self.search_patterns):
                txt = parse.flatten()
                found_obj = self.search_proc(txt, pattern, 0)
                rule_identifier = self.name + ('_pat_' + str(pattern_num) if len(self.search_patterns) > 1 else '')
                while found_obj:
                    start = found_obj.start()
                    end = found_obj.end()
                    left_context = txt[:start]
                    right_context = txt[end:]
                    new_strucs = self.replacement_proc(found_obj, left_context, right_context, parse)
                    if new_strucs is None or end <= start:
                        # Replacement_proc did not validate what the search_proc found.
                        end += 1
                    elif self.change_parse_in_place:
                        parse.resegment(new_strucs = new_strucs, start_pos_new_strucs = start, end_pos_new_strucs = end - 1, rule_name = rule_identifier)
                        parse_changed_or_new_parse_found = True
                        parse.rules_utilized.append(rule_identifier)
                        parse.changed_on_last_pass = True
                        end = 0
                    else:
                        # use copy, not deepcopy, because we specifically want to preserve references to strucs.
                        #new_parse = copy.copy(parse)
                        #new_parse.strucs = copy.copy(parse.strucs)
                        new_parse = parse.copy()
                        new_parse.resegment(new_strucs = new_strucs, start_pos_new_strucs = start, end_pos_new_strucs = end - 1, rule_name = rule_identifier)

                        is_new_parse = True
                        # Check if the new parse is pointwise different from other parses
                        for a_parse in parse.sig.parses + new_parses:
                            if new_parse.pointwise_equal_by_strucs(a_parse):
                                is_new_parse = False
                                break
                        if is_new_parse:
                            new_parse.rules_utilized.append(rule_identifier)
                            new_parses.append(new_parse)
                            parse_changed_or_new_parse_found = True
                            new_parse.changed_on_last_pass = True


                    txt = parse.flatten()
                    found_obj = self.search_proc(txt, pattern, end)     # specify to search after end char because if we are not changing parse in place
                                                                        # the actual parse.flatten() remains unchanged throughout the iterations of the loop.

            # Delete redundant parses, i.e parses where the parse.strucs lists are pointwise identical.
            if parse_changed_or_new_parse_found and self.change_parse_in_place:
                for a_parse_num, a_parse in enumerate(parse.sig.parses):
                    if a_parse_num != parse_num and parse.pointwise_equal_by_strucs(a_parse):
                        del parse.sig.parses[parse_num: parse_num + 1]
                        break

        sig.parses.extend(new_parses)



##### UTILITIES #######
def trim(text):
    """    Trims white spaces  at the front, end, and double spaces inside a string    """

    pattern_trim_1 = re.compile('^\s+')         # at beginning: remove leading spaces
    pattern_trim_2 = re.compile('\s+$')         # at end: remove trailing spaces
    pattern_trim_3 = re.compile('\s+\.$')       # at end of string: xxx . --> xxx.
    pattern_trim_4 = re.compile('\s\s+')        # remove double space

    text = text.strip()
    text = pattern_trim_1.sub('', text)
    text = pattern_trim_2.sub('', text)
    text = pattern_trim_3.sub('.', text)
    text = pattern_trim_4.sub(' ', text)

    return text

def normalize_string(text):
    """Initial normalizaiton of raw sig: lowercases and does basic character cleanup. """

    text = text.lower()
    text = re.sub(r'\n', ' ', text)              # replace carriage returns with spaces
    text = re.sub(r'_?x000d_?', ' ', text)       # replace hex carriage returns with spaces
    text = re.sub(r'(\[|{|\()\s*', '(', text)    # replace brackets with parentheses
    text = re.sub(r'\s*(\]|}|\))', ')', text)    # replace brackets with parentheses
    text = re.sub(r'\*', ' ', text)              # remove stars
    text = re.sub(r'\|', ' ', text)              # remove vertical bars
    text = re.sub(r'"', ' ', text)
    text = re.sub(r"`", "'", text)
    text = re.sub(r'&', ' and ', text)           # & --> "and"
    text = re.sub(r'!+', '.', text)              # ! -> .
    text = re.sub(r'\s+,', ',', text)            # " ," -> ","
    text = re.sub(r',(?=[a-z])', ', ', text)     # ",xyz" -> ", xyz"
    text = re.sub(r';', ' ', text)               # ";" -> " " (";" is often used in Latin sig abbrs)
    #text = re.sub(r'(?<=[a-z0-9.,])\(', ' (', text)  # "x(" -> " ("  (run on open paren)
    #text = re.sub(r'\)(?=[a-z0-9])', ') ', text)  # ")x" -> ") "  (run on close paren)
    text = re.sub(r'(?<![a-z])once(?![a-z])', r'1 time', text)        # once -> 1 time
    text = re.sub(r'(?<![a-z])twice(?![a-z])', r'2 times', text)      # twice -> 2 times
    text = re.sub(r'(?<![a-z])thrice(?![a-z])', r'3 times', text)     # thrice -> 3 times
    text = re.sub(r'(?<![0-9])1st(?![a-z])', r'first', text)          # 1st -> first
    text = re.sub(r'(?<![0-9])2nd(?![a-z])', r'second', text)         # 2nd -> second
    text = re.sub(r'(?<![0-9])3rd(?![a-z])', r'third', text)          # 3rd -> third
    text = re.sub(r'(?<![a-z])(each day|everyday)(?![a-z])', r'every day', text)         # each day -> every day
    text = re.sub(r'(?<![a-z])nightly(?![a-z])', r'every night', text)        # nightly -> every night
    text = re.sub(r'(?<![a-z])on day (one|1) and(?![a-z])', r'on day 1, and', text)   # inserts a comma to prevent incorrect parsing of phrases such as
                                                                                    # "take 1 tspoon on day one and 1/2 teaspoon on days 2-5" by seeing "one and 1/2" as 1.5,
                                                                                    # as opposed to "take 1 tspoon on day 1, and take 1/2 tspoon on days 2-5"
    # remove duplicates
    text = re.sub(r'(?<![a-z])apply apply(?![a-z])', r'apply', text)
    text = re.sub(r'(?<![a-z])instill? instill(?![a-z])', r'instill', text)
    text = re.sub(r'(?<![a-z])inject inject(?![a-z])', r'inject', text)
    text = re.sub(r'(?<![a-z])insert insert(?![a-z])', r'insert', text)
    text = re.sub(r'(?<![a-z])take take(?![a-z])', r'take', text)
    text = re.sub(r'(?<![a-z])use use(?![a-z])', r'use', text)
    text = re.sub(r'(?<![a-z])tablet tablet(?![a-z])', r'tablet', text)
    text = re.sub(r'(?<![a-z])teaspoonful teaspoonful(?![a-z])', r'teaspoonful', text)
    text = re.sub(r'(?<![a-z])for for(?![a-z])', 'for', text)
    text = re.sub(r'(?<![a-z])every every(?![a-z])', 'every', text)
    text = re.sub(r'(?<![a-z])in each nostrill?s? (in)?(to)?\s?each nostrill?s?(?![a-z])', r'in each nostril', text)    # frequent duplicaiton of "each nostril"
    text = re.sub(r'(?<![a-z])each nostrill?s? daily each nostrill?s? daily(?![a-z])', r'each nostril daily', text)    # frequent duplicaiton of "each nostril"

    text = re.sub(r'(?<![a-z])\(1 fl oz \= 30 ml\)(?![a-z])', r' ', text)    # remove from dictionary an explanatory note "(1 FL OZ = 30 ML)"


    text = trim(text)   # clean up white space
    return text

def spelling_corrections(text):
    """ Correction of common misspellings that don't need context for dismbiguation (to the best of the data we have seen so far)
    """

    text = re.sub(r'(?<![a-z])affeceted|afected|afective|afeected(?![a-z])', r'affected', text)
    text = re.sub(r'(?<![a-z])aply|applied|appy|aplly|aaply|aapply|apppy(?![a-z])', r'apply', text)
    text = re.sub(r'(?<![a-z])ashtma|astma|asmatha(?![a-z])', r'asthma', text)
    text = re.sub(r'direceted', r'directed', text)
    text = re.sub(r'daily for (\d+) daily', r' daily for \1 days ', text)
    text = re.sub(r'\s?daily\s?', r' daily ', text)
    text = re.sub(r'docotr|docotor', r'doctor', text)
    text = re.sub(r'(?<![a-z])(evry|ever)(?![a-z])', r'every', text)
    text = re.sub(r'effected', r'affected', text)
    text = re.sub(r'(?<![a-z])eyes? lid', r'eyelid', text)        # "eye lids" --> "eyelid" mostly to help us differentiate
    text = re.sub(r'f00d', r'food', text)
    text = re.sub(r'(?<![a-z])inahle', r'inhale', text)
    text = re.sub(r'(?<![a-z])instil(?![a-z])', r'instill', text)
    text = re.sub(r'instrution|instructuion', r'instruction', text)
    text = re.sub(r'(tablets?|capsules?)\s(my|buy|mu|bu|per|ny) mouth', r'\1 by mouth', text)
    text = re.sub(r'more then', r'more than', text)
    text = re.sub(r'physici?an', r'doctor', text)
    text = re.sub(r'(?<![a-z])(nabulizer|nebuliser|nebuliver|nebulzer|nebulier|nebuliazer)(?![a-z])', r'nebulizer', text)
    text = re.sub(r'to relieve\s?pain', r'for pain', text) # because of the frequent misspelling of "relive" for "relief" (if no relief call 911) need to get rid of frequent legitimate uses
    text = re.sub(r'to relieve\s?irritation', r'to reduce irritation', text) # because of the frequent misspelling of "relive" for "relief" (if no relief call 911) need to get rid of frequent legitimate uses
    text = re.sub(r'(?<![a-z])relive?|relife?|relieve?(?![a-z])', r' relief ', text)        # "relief" or "relieve" misspelling (but even if "relieve" we utilize it as "if no relief ..."
    text = re.sub(r'suppliment', r'supplement', text)
    text = re.sub(r'supposotory', r'suppository', text)
    text = re.sub(r'(?<![a-z])(spary|sprae|sprasy|sprayd|spraye|sprray|spyary|sray|srpay)(s?)(?![a-z])', r'spray\2', text)
    text = re.sub(r'tabelt', r'tablet', text)


    text = re.sub(r'(?<=call)\s?9\s?\-\s?1\s?\-\s?1\s?(?![0-9])', r' 911 ', text)    # replace "call 9-1 - 1" with "call 911"
    text = re.sub(r'(?<![a-z])dialy(?![a-z])', r'daily', text)
    text = re.sub(r'(?<![a-z])dayly(?![a-z])', r'daily', text)
    text = re.sub(r'(?<![a-z])aday(?![a-z])', r'a day', text)
    text = re.sub(r'(?<![a-z])s day(?![a-z])', r'a day', text)
    text = re.sub(r'(?<=\d)daily(?![a-z])', r' daily', text)                         #missing space before "daily", e.g. "1daily"
    text = re.sub(r'(?<![a-z])dialyis|dialysys(?![a-z])', r'dialysis', text)


    text = re.sub(r'(\d+)\s?g(?![a-z])', r'\1 gram ', text)                         # 17g means 17 gram if not followed by other letters
    text = re.sub(r'>>\sgm?(?=(\s|\.))', r'>> gram', text)                          # For matching dictionary entries for Gram.
    text = re.sub(r'(\d+)l(?![a-z])', r'\1 liter ', text)                           # 2L means 2 liters if not followed by other letters

    text = re.sub(r'(?<![a-z])see (accompanying )?insert(?![a-z])', r' see attached sheet ', text)      # we need to distinguish "insert" as directive vs "insert" as noun.
    text = re.sub(r'(?<![a-z])(package|product|attached|accompanying|with|on|an) insert(?![a-z])', r' attached sheet ', text)      # we need to distinguish "insert" as directive vs "insert" as noun.


    text = trim(text)   # clean up white space
    return text

def latin_2_sig(text):
    """ Expands Latin and pseudo latin abbreviations within the sig.

        Returns expanded_sig.
        Brian's standard examplesa are 'ap bid F10d UD' and 'T1TPOTIDPRN'
    """

    def substitution_loop(string, loop_num):
        """ We need to do a recurive loop because some of disabbreviations are sensitive to surrounding spaces, while
            others are absolute. So by doing one round of substitutions you may enable new substitutions to be made
            in the next round.
        """
        source_string = string

        # Form, Directive, and Dose quantity
        string = re.sub(r'(?<![a-z])t(\d+)tsp(?![a-z])', r' take \1 teaspoon ', string)
        string = re.sub(r'(?<![a-z])t(\d+)tb(sp?)?(?![a-z])', r' take \1 tablespoon ', string)
        string = re.sub(r't1t\.?', r' take 1 tablet ', string)
        string = re.sub(r't(\d\.?\d?)t', r' take \1 tablet ', string)
        string = re.sub(r'(\d+)tpo', r' \1 tablets by mouth ', string)

        string = re.sub(r't1c\.?', r' take 1 capsule ', string)
        string = re.sub(r't(\d\.?\d?)c', r' take \1 capsule(s) ', string)
        string = re.sub(r'1cpo', r' 1 capsule by mouth ', string)
        string = re.sub(r'(\d+)cpo', r' \1 capsules by mouth ', string)

        string = re.sub(r'(^|\s)(\d+)(capsules?|c(?![a-z]))', r' \2 capsule ', string)
        string = re.sub(r'(\d+) c(?![a-z])', r' \1 capsule ', string)
        string = re.sub(r'(^|\s)(\d+)(tablets?|t(?![a-z]))', r' \2 tablet ', string)
        string = re.sub(r'(\d+) t(?![a-z])', r' \1 tablet ', string)

        string = re.sub(r'^(\d)(?=q)', r' \1 tablet ', string)              # "2q.." at the beginning of sig means "2 tablets every..."
        string = re.sub(r'^(\d)\-(\d)(?=q)', r' \1-\2 tablets ', string)    # "1-2q.." at the beginning of sig means "1-2 tablets every..."

        string = re.sub(r'gtts', r' drops ', string)
        string = re.sub(r'gtt', r' drop ', string)
        string = re.sub(r'(?<![a-z])supp?(?![a-z])', r' suppository ', string)
        string = re.sub(r'(\d+)\s?u(?![a-z])', r' \1 units ', string)

        string = re.sub(r'^tk?\s', r'take ', string)
        string = re.sub(r'^tk?(\d)', r'take \1', string)
        string = re.sub(r'^ap?\s', r'apply ', string)
        string = re.sub(r'^ap?(\d)', r'apply \1', string)
        string = re.sub(r'^u\s', r'use ', string)
        string = re.sub(r'^u(\d)', r'use \1', string)
        string = re.sub(r'^il(?![a-z])', r' instill ', string)
        string = re.sub(r'(^|\s)u as directed', r' use as directed', string)
        string = re.sub(r'(^|\s)aaa(?![a-z])', r' apply to affected area ', string)

        # x and q (Duration & Freq)

        # xN generally means "for N (days, etc)" but can mean "times", e.g. "repeat x2" = repeat twice.
        # Nx generally means "N times"
        string = re.sub(r'(x|f)(\d+)d(ays?|ys)?', r' for \2 days ', string)
        string = re.sub(r'(x|f)(\d+)w(eeks?|ks?)?', r' for \2 weeks ', string)
        string = re.sub(r'(x|f)(\d+)h(ours?|rs?)?', r' for \2 hours ', string)
        string = re.sub(r'(?<![a-z])x\s?(\d+)\s?d(ays?|ys)?', r' for \1 days ', string)         # "x 7 days"
        string = re.sub(r'(?<![a-z])x\s?(\d+)\s?(\-|to)\s?(\d+)\s?d(ays?|ys)?', r' for \1 to \3 days ', string)    # "x 3-5 days"
        string = re.sub(r'(?<![a-z])x\s?(\d+)\s?w(eeks?|ks?)?', r' for \1 weeks ', string)      # "x 7 weeks"
        string = re.sub(r'(?<![a-z])x\s?(\d+)\s?$', r' for \1 ', string)                        # at the end of the sig, "x7" means "every 7 timeunits" where timeunits are implied before.
        string = re.sub(r'repeat x(\d+)', r' repeat \1 times ', string)
        string = re.sub(r'(?<!every)(?<!q)\s(\d+)d(\s|$)', r' for \1 days ', string)                  # " 7d " means "for 7 days" unless preceded by "every'
        string = re.sub(r'(\d+)x', r' \1 times ', string)                                       # 7x daily
        string = re.sub(r'(\d+) x(\s|$)', r' \1 times ', string)                                # 7x daily

        string = re.sub(r'bid', r' 2 times a day ', string)
        string = re.sub(r'tid', r' 3 times a day ', string)
        string = re.sub(r'qid', r' 4 times a day ', string)

        string = re.sub(r'qam|(^|\s)q am|every a.?m.?', r' every morning ', string)
        string = re.sub(r'qpm|(^|\s)q pm|every p.?m.?', r' every evening ', string)

        # q means "every" (except in qid). q4h = every 4 hours, qd = every day, q8-12h = every 8-12 hours

        # q hours.
        # Note: "q72" means every 72 hours but "q46" and "q46h" means "every 4-6 hours" even if not followed by "h". I.e. dash can be omitted depending on the numbers.
        # Note: "q4" means "every 4 HOURS" even if "h" is not present. Same for q3, q4, q6, q8, q24 and q72.


        string = re.sub(r'q(\d+)(to|\-)(\d+)h(ours?|rs?)?(\.|\s|$)', r' every \1 to \3 hours ', string)
        string = re.sub(r'(^|(?<=\d|\s))q46(h(ours?|rs?)?)?(\.|\s|$)', r' every 4 to 6 hours ', string)     # q46 before space = every 4-6 hours even without "h" or dash
        string = re.sub(r'(^|(?<=\d|\s))q68(h(ours?|rs?)?)?(\.|\s|$)', r' every 6 to 8 hours ', string)     # q68 before space = every 6-8 hours even without "h" or dash
        string = re.sub(r'q46h(ours?|rs?)?', r' every 4 to 6 hours ', string)                 # q46h = every 4-6 hours even without dash
        string = re.sub(r'q68h(ours?|rs?)?', r' every 6 to 8 hours ', string)                 # q68h = every 6-8 hours even without dash
        string = re.sub(r'(^|(?<=\d|\s))q(2|3|4|6|8|24|72)(h(ours?|rs?)?)?(\.|\s|$)', r' every \2 hours ', string) # q72 = every 72 hours even without an "h". Same q24, q4, q6
        string = re.sub(r'(^|(?<=\d|\s))q12\s?$', r' every 12 hours ', string)                      # q12 terminal = every 12 hours
        string = re.sub(r'q(\d+)h(ours?|rs?)?', r' every \1 hours ', string)      # q72h = every 72 hours
        string = re.sub(r'q(\d+) (hours?|hrs?)(\.|\s|$)', r' every \1 hours ', string)

        string = re.sub(r'(\d+)\s?h(rs?)?(,|\s|$)', r' \1 hours ', string)

        # q minutes
        string = re.sub(r'q(\d+)\s?m(ins?)?(\.|\s|$)', r' every \1 minutes ', string)

        # q days. qod = every other day. q2d = every 2 days
        string = re.sub(r'qd|(^|\s)q d(ays?)?(\.|\s|$)', r' every day ', string)
        string = re.sub(r'q(\d+)d(ays?)?', r' every \1 days ', string)
        string = re.sub(r'(^|(?<=\d|\s))qod(\s|$)', r' every 2 days ', string)

        string = re.sub(r'(\d+)\s?d(,|\s|$)', r' \1 days ', string)


        # q weeks
        string = re.sub(r'qwk|(^|\s)q (wk|week)s?(\.|\s|$)', r' every week ', string)
        string = re.sub(r'qw(\s|$)', r' every week ', string)
        string = re.sub(r'q(\d+)w(eeks?|ks?)?', r' every \1 weeks ', string)  # q2h = every 2 weeks

        #catch all q
        string = re.sub(r' q ', r' every ', string)




        # Route
        string = re.sub(r'(^|(?<=(\d|\s)))po(,|\s|$)', r' by mouth ', string)
        string = re.sub(r'po\s?(?=every|q(\d|d|am|pm)|\d+ times|daily)', r' by mouth ', string)    # even if not preceded by space, allow if followed by q\d
        string = re.sub(r'(^|(?<=\d|\s))(sc|sq)(,|\.)?(\s|$)', r' subcutaneously ', string)
        string = re.sub(r'(^|(?<=\d|\s))(im)(,|\.)?(\s|$)', r' intramuscularly ', string)
        string = re.sub(r'(^|(?<=\d|\s))pv(,|\.)?(\s|$)', r' vaginally ', string)

        # Timing
        string = re.sub(r'(^|\s)pc(\s|$)', r' after meals ', string)
        string = re.sub(r'(^|\s)ac(\s|$)', r' before meals ', string)
        string = re.sub(r'(^|\s)wf(\s|$)', r' with food ', string)
        string = re.sub(r'(^|\s)wc|wm(\s|$)', r' with meals ', string)

        string = re.sub(r'(\d)hs', r' \1 at bedtime ', string)
        string = re.sub(r'(\d)ths', r' \1 tablet at bedtime ', string)
        string = re.sub(r'(\d)qhs', r' \1 every day at bedtime ', string)
        string = re.sub(r'qhs|(^|\s)q hs', r' every day at bedtime ', string)
        string = re.sub(r'(^|(?<=\d|\s))hs(\s|$)', r' at bedtime ', string)

        string = re.sub(r'q(10|11|12|\d)pm(?![a-z])', r' every day at \1pm ', string)    # q8pm means "every day at 8 pm"
        string = re.sub(r'q(10|11|\d)am(?![a-z])', r' every day at \1am ', string)

        string = re.sub(r'(?<![a-z0-9])in (the )?a\.?m.?(?![a-z])', r' in the morning ', string)
        string = re.sub(r'(?<![a-z0-9])(?<!\d\s)a\.?m.?(?![a-z])', r' in the morning ', string)
        string = re.sub(r'(?<![a-z0-9])in (the )?p\.?m.?(?![a-z])', r' in the evening ', string)
        string = re.sub(r'(?<![a-z0-9])(?<!\d\s)p\.?m.?(?![a-z])', r' in the evening ', string)


        # Site
        string = re.sub(r'(^|\s)ou(\s|$)', r' each eye ', string)
        string = re.sub(r'(^|\s)os(\s|$)', r' left eye ', string)
        string = re.sub(r'(^|\s)od(\s|$)', r' right eye ', string)
        string = re.sub(r'(\d+)gou(\s|$)', r' \1 drop in each eye ', string)

        # Misc
        string = re.sub(r'prn', r' as needed ', string)
        string = re.sub(r' p r n(?![a-z])', r' as needed ', string)
        string = re.sub(r'(?<![a-z])u?ud|utd(?![a-z])', r' use as directed ' , string)
        string = re.sub(r'(^|\s)ug(\s|$)', r' until gone ', string)
        string = re.sub(r'(^|\s)uf(\s|$)', r' until finished ', string)
        string = re.sub(r'(?<![a-z])stat(?![a-z])', r' now ', string)

        expanded_string = trim(string)

        if debug and loop_num > 5:
            raise Exception('Finished %d loops in latin_2_sig substitution_loop(). Original Source: ||%s|| Tranformed to: ||%s||' % (loop_num, original_string, expanded_string))
        elif debug:
            # print('Finished loop %d in latin_2_sig substitution_loop(). Source: ||%s|| Tranformed to: ||%s||' % (loop_num, source_string, expanded_string))
            pass
        if expanded_string == source_string:
            return expanded_string
        else:
            loop_num += 1
            return substitution_loop(expanded_string, loop_num)

    original_string = text
    return substitution_loop(text, loop_num = 0)



def trim_dictionary(a_dict):
    """Trims all items in the dictionary that have None or '' values. """

    return dict([(key, val) for (key, val) in a_dict.items() if val])

def padlines(txt, width):
    """ Pads each line of txt on the left with spaces of width width."""

    return '\n'.join([' ' * width + line for line in txt.splitlines()])


def get_struc_labels(struc_list, delimiter = None, omit_spaces_and_punctuation = None):
    delimiter = '' if delimiter is None else delimiter

    if omit_spaces_and_punctuation:
        struc_list = [struc for struc in struc_list if not struc.is_space_or_punctuation_only()]

    struc_list_string = delimiter.join([struc.label for struc in struc_list])
    return struc_list_string

def reset_list(alist):
    """ Resets list to empty list without reassignment. Useful for resetting global variable errors_list without declaring it global
    """

    del alist[0:len(alist) + 1]


def cross_lists(*lists):
    """ Takes several lists and returns a modified cross product: a list of all lists of type [x0, x1, ..., xn] where xi is an element from xi'th list.
        E.g., if list1 = [1, 2] and list2 = ['a', 'b'], cross_lists(list1, list2) = [[1, 'a'], [1, 'b'], [2, 'a'], [2, 'b']]
    """

    def cross2(list1, list2):
        rows = [el1 + [el2] for el1 in list1 for el2 in list2]
        return rows

    if len(lists) == 1:
        return [[el] for el in lists[0]]
    last_list = lists[-1]
    initial_lists = lists[:-1]
    cross_initial_lists = cross_lists(*initial_lists)
    final_cross = cross2(cross_initial_lists, last_list)
    return final_cross

def choose_m_of_n(alist, m):
    """ From a list alist choose all sublists of length m. Return the list of these sublists.

        Each sublist remains in the same order as the original alist.

        Example: alist = [1, 2, 3], m = 2. Returns [[1, 2], [2,3], [1,3]]
    """

    if m < 1 or m > len(alist):
        return []
    elif len(alist) == m:
        return [alist]
    elif m == 1:
        return [[struc] for struc in alist]
    else:
        element = alist[-1]
        rest = alist[:-1]
        without_element = choose_m_of_n(rest, m)
        with_element = [sublist + [element] for sublist in choose_m_of_n(rest, m-1)]
        result = with_element + without_element
        return result



def instantiate_atom_string_numerically(atom_string, numeric_mappings, locale):

    for (variable, value) in numeric_mappings:
        if locale not in locales_using_period_for_decimal:
            value = value.replace('.', ',')
        atom_string = atom_string.replace(variable, value)
    return atom_string

##### END UTILITIES #######


def rule_numerics_identification():

    def search_proc(txt, search_pattern = None, start = None):
        """ search_patterns is a dummy parameter used for consistent interface with rule.apply_rule()
        """

        numerics_segment = numeric_id.extract_numerical_expressions(txt, start = start, return_first_found = True)
        return numerics_segment

    def replacement_proc(num_seg, left_context, right_context, parse):

        seg_type = num_seg.label

        if seg_type == 'DATE':
            quant_struc = Date(day = num_seg.day, month = num_seg.month, year = num_seg.year)
        elif seg_type == 'MAXDOSE':
            quant_struc1 = Quant(num_type = 'int', value = num_seg.dose)
            quant_struc1.constituents = [str(num_seg.dose)]
            quant_struc2 = Quant(num_type = 'int', value = num_seg.hours)
            quant_struc2.constituents = [str(num_seg.hours)]
            return [quant_struc1, Struc(label = '/', constituents = '/'), quant_struc2]
        elif seg_type == 'BLOODPRESSURE':
            quant_struc = BloodPressureRatio(systolic = num_seg.systolic, diastolic = num_seg.diastolic)
        elif seg_type == 'FLOAT' and num_seg.float_type == 'FRAC':
            quant_struc = FracQuant(integer = num_seg.integer, numerator = num_seg.numerator, denominator = num_seg.denominator)
        else:
            quant_struc = Quant(num_type = seg_type.lower(), value = num_seg.value)
        quant_struc.constituents = [num_seg.replaced_txt]
        return [quant_struc]

    rule = Rule_ExtractStrucs(   name = 'numerics_identification',
                                    search_patterns = [''],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['preprocessing']),
                                    change_parse_in_place = True
                                    )
    return rule



def rule_assembly_form():

    pattern = re.compile(r'''
                                    (?<![a-z])                                  # not preceded by a letter
                                    (the \s)?
                                    (                                           # Form type: one of the following:
                                        (?P<application>application)    |       # application
                                        (?P<applicatorful>applicator(fool|full?))       |   # applicatorful
                                        (?P<bottle>bottle)              |       # bottle
                                        (?P<can>can)                    |       # can (to drink contents of)
                                        (?P<capful>cap(fool|full?))     |       # capful
                                        (?P<cc>cc)                      |       # CC
                                        (?P<capsule>
                                            gel \s? cap(sule)?  |
                                            cap(sule)?
                                        )                               |       # capsule
                                        (?P<dropperful>
                                            dropp?er(fool|full?)
                                        )                               |       # dropperful
                                        (?P<drop>drop)                  |       # drops
                                        (?P<gram>gram|gm)               |       # grams
                                        (?P<inhalation>                         # inhalation
                                            (?<!for \s)                         # "for inhalation" is really directive "inhale"
                                            inhalation
                                        )                               |
                                        (?P<liter>liter)                |        # liter
                                        (?P<lozenge>lozenge)            |       # lozenges
                                        (?P<mcg>mcg)                    |       # MCG
                                        (?P<mg>mg)                      |       # MG
                                        (?P<ml>                                 # ML
                                            ml \s? \( cc \)     |               # "ml(cc)"  which we take to mean just "ml"
                                            milliliter          |
                                            ml
                                         )                              |       #
                                        (?P<ounce>oz|ounce)             |       # Ounces
                                        (?P<packaging>packaging)        |       # packaging
                                        (?P<packet>
                                            packet      |
                                            pkt
                                        )                               |       # packet
                                        (?P<pack>pack)                  |       # pack (seems to have multiple meanings, default to directive "use")
                                        (?P<pad>pad)                    |       # pad
                                        (?P<patch>patch)                |       # patch
                                        (?P<puff>puff)                  |       # puff
                                        (?P<ring>ring)                  |       # vaginal ring
                                        (?P<scoop>                              # scoop
                                            scoop(full?)?
                                        )                               |
                                        (?P<spray>spray)                |       # spray
                                        (?P<suppository>
                                            supp?ositor(y|ies)
                                        )                               |       # suppository
                                        (?P<syringe>syringe)            |       # syringe
                                        (?P<tablet>                             # tablet (including pill, as long as it is not "water pill", which is synonym for diuretic and is handled in rule_miscellanous_diuretic()
                                            tab(let)?           |
                                            (?<!water\s) pill
                                        )                               |
                                        (?P<teaspoon>tea?spoon(s?fool|s?full?)?|tsp)         |   # teaspoon
                                        (?P<tablespoon>table \s? spoon(s?fool|s?full?)?|tbsp)    |   # tablespoon
                                        (?P<vial>
                                            vial        |
                                            vail        |
                                            respule     |
                                            ampoule?    |
                                            ampull?e    |
                                            amp
                                        )                               |     # vial
                                        (?P<unit>                               # unit
                                            unit
                                        )
                                    )
                                    (                                           # possibly followed by indication of plurality
                                        (?P<plural>e?s | \'s)   |                     # simple s or es for plurality
                                        (?P<plurality_either>                   # or indication of unspecified plurality such as tablet/s or tab(s)
                                            \/s         |
                                            \(e?s\)
                                        )
                                    )?
                                    \.?
                                    (?![a-z])                                   # not followed by another letter
                                ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        groupname2group = trim_dictionary(match_object.groupdict())

        if 'plural' in groupname2group:
            plurality = 'plural'
        elif 'plurality_either' in groupname2group:
            plurality = 'plurality_either'
        else:
            plurality = 'singular'
        # the form type is that group name found that is not a group name for plurality identifier group.
        form_types = set(groupname2group.keys()) - Form.plurality_permissible_values
        if len(form_types) > 1 and debug:
            raise Exception('More than 1 pill form found in word %s' % txt)
        form_name = list(form_types)[0]

        if form_name in ('mcg', 'mg', 'ml', 'ounce', 'liter'):        # make sure that short abbreviations such as mg are preceded by a number
            left_context_pattern = re.compile('QUANT\s?\-?\s?$')
            if not left_context_pattern.search(left_context):
                return None

        if form_name == 'spray' and len(left_context) < 2:      # "spray" can be a verb, and if it occurs at the beginning of the string, it IS a directive
            return None

        if form_name == 'pack':
            if 'mix' in left_context or 'dissolve' in left_context or 'drink' in left_context:
                form_name = 'packet'
            elif 'DIRECTIVE' in left_context:
                left_context_pattern = re.compile('(?P<directive>DIRECTIVE)')
                left_directive_obj = left_context_pattern.search(left_context)
                if left_directive_obj:
                    left_directive_start = left_directive_obj.start('directive')
                    left_directive = parse.position2struc(left_directive_start)
                    if left_directive.value in ('mix', 'dissolve', 'drink'):
                        form_name = 'packet'
            else:
                # no directive. E.g. a frequent full sig: "1 pack once daily'
                form_name = 'packet'
            if form_name == 'packet':
                # we changed the form from pack to packet. Record this rules_used.
                form = Form(form_name = form_name, plurality = plurality, constituents = [match_object.group()])
                form.rules_used.append('*deduced_from_pack_in_form_identification*')


        form = Form(form_name = form_name, plurality = plurality, constituents = [match_object.group()])

        return [form]

    rule = Rule_ExtractStrucs(   name = 'form',
                                    search_patterns = [pattern],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule





def rule_assembly_directive():

    pattern = re.compile(r'''
                            (?<![a-z])
                            (use \s (to|for) \s)?                 # frequent prefix "use to test blood sugar" or "use to inject". Merge with the directive
                            (for \s)?                       # e.g. "for testing", "for use with"
                            (?P<directive>
                                apply                       |
                                (?P<check>
                                    check(ing)?
                                )            |
                                chew \s (and|AND_CONJ) \s swallow      |
                                chew                        |
                                (?P<drink>
                                    ((and|AND_CONJ|then|THEN_CHRONO) \s?)?
                                    drink
                                    (
                                        (\s the)?
                                        \s contents?
                                        (\s of)?
                                    )?
                                )                           |
                                (?P<dissolve>                   # "dissolve & take" is dissolve.
                                    diss?olve
                                    (\s (and|AND_CONJ|then|THEN_CHRONO) \s take)?
                                )                           |
                                give                        |
                                (?P<inhale>                     # for some reason expressions such as "take 2 ml 2 times a day inhaled" where "inhaled" is really a directive to replace "take" and "use"
                                    inhaled?            |
                                    for \s inhalation           # eg "use 1 capsule for inhalation". We transform "for inhalation" into "inhale" and later change "use" to "inhale"
                                )                           |
                                inject                      |
                                (?P<insert>
                                    (?<!unwrap\sand\s)
                                    insert
                                )                           |
                                instill                     |
                                (?P<mix>
                                    mix (\s (and|AND_CONJ) \s drink)    |
                                    mix
                                )                           |
                                (?P<place>
                                    (?<!in\s)                   # avoid "leave in place" = "remove"
                                    place
                                )                           |
                                put                         |
                                (?P<rinse>
                                            rinse
                                            (\s out)?
                                            (\s off)?
                                            (\s with)?
                                            ((and|AND_CONJ|then|THEN_CHRONO) \s spit (\s out)? )?
                                )                           |
                                shampoo                     |
                                spray                       |
                                (?P<stop>                       # used in "take for 3 weeks then stop (for 1 week)?"
                                    stop            |
                                    discontinue
                                )                           |
                                (?P<swallow>
                                    (?<!swish \s and \s)
                                    (?<!swish \s AND_CONJ \s)
                                    swallow
                                )                           |
                                (?P<swish_and_spit>
                                    swish \s (and|AND_CONJ|then|THEN_CHRONO) \s spite?
                                )                           |
                                (?P<swish_and_swallow>
                                    swish \s (and|AND_CONJ) \s swallow
                                )                           |
                                (?P<swish>
                                    swish
                                )                           |
                                (?P<take>take (?!\s off) )  |   # "take off" = "remove"
                                (?P<test>
                                    (
                                        (use \s for \s)? testing   |
                                        (to \s)? test
                                    )
                                )                           |
                                (?P<unwrap_insert>
                                    unwrap
                                    \s?
                                    (and|AND_CONJ)?
                                    \s?
                                    insert
                                )                           |
                                unwrap                      |
                                use                         |
                                (?P<remove>                                 # remove is last pattern because it can capture "off", which is also part of "rinse off"
                                    (remove \s (and|AND_CONJ) \s)? leave \s off     |
                                    (remove \s (and|AND_CONJ) \s)? leave \s out     |
                                    remove              |
                                    peel \s off         |
                                    take \s off         |
                                    (?P<off_alone>off)
                                )
                            )
                            (?P<test_sugar>                 # optional subject of test: "test blood sugar", "check glucose level", etc.
                                \s
                                (for \s)?
                                (blood \s)?
                                (sugar|glucose)      # mandatory for test_sugar
                                (\s levels?)?
                            )?
                           (?![a-z])
                           ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())

        directive = groupname2group['directive']
        dir_attribute = None

        if 'check' in groupname2group:
            directive = 'test'
        elif 'dissolve' in groupname2group:
            directive = 'dissolve'
        elif 'drink' in groupname2group:
            # There are 2 types of cases with "drink". Most typically, it is used as secondary directive after "mix" or "dissolve": "mix 17 grams in 8 oz liquid and drink twice daily"
            # Here "drink" is not a separate directive that requires a separate Event, but an instantaneous follow-up to the main event. We have to get rid of it.
            # But a bona-fide use of "drink" is "drink 2 cans 3 times daily". These we have to keep.
            directive = 'drink'
            left_context_pattern = re.compile('(?P<directive>DIRECTIVE)')
            left_directive_obj = left_context_pattern.search(left_context)
            if left_directive_obj:
                return []
        elif 'inhale' in groupname2group:
            # for some reason expressions such as "take 2 ml 2 times a day inhaled" where "inhaled" is really a directive to replace "take" and "use"
            directive = 'inhale'
            left_context_pattern = re.compile('(?P<directive>DIRECTIVE)')
            left_directive_obj = left_context_pattern.search(left_context)
            if left_directive_obj:
                left_directive_start = left_directive_obj.start('directive')
                left_directive = parse.position2struc(left_directive_start)
                if left_directive.value in ('take', 'use'):
                    left_directive.value = 'inhale'
                    return []
        elif 'insert' in groupname2group:
            directive = 'insert'
        elif 'mix' in groupname2group:
            directive = 'mix'
        elif 'place' in groupname2group:
            directive = 'place'
        elif 'remove' in groupname2group:
            directive = 'remove'
        elif 'rinse' in groupname2group:
            directive = 'rinse'
        elif 'stop' in groupname2group:
            directive = 'stop'
        elif 'swallow' in groupname2group:
            directive = 'swallow'
        elif 'swish' in groupname2group:
            directive = 'swish'
        elif 'swish_and_swallow' in groupname2group:
            directive = 'swish and swallow'
        elif 'swish_and_spit' in groupname2group:
            directive = 'swish and spit'
        elif 'take' in groupname2group:
            directive = 'take'
        elif 'test' in groupname2group:
            directive = 'test'
        elif 'unwrap_insert' in groupname2group:
            directive = 'unwrap and insert'

        if 'test_sugar' in groupname2group:
            dir_attribute = 'blood_sugar'
            if directive == 'check':
                directive = 'test'
        dir_struc = Directive(value = directive)
        dir_struc.constituents = [match_object.group()]

        if dir_attribute:
            dir_attr_struc = DirAttribute(value = dir_attribute)
            space = Struc(label = ' ')
            return [dir_struc, space, dir_attr_struc]
        else:
            return [dir_struc]

    rule = Rule_ExtractStrucs(   name = 'directive',
                                    search_patterns = [pattern],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule

def rule_assembly_broken_up_directives():
    pattern0 = re.compile(r'''      # Broken complex directive. E.g. "Chew 1 tablet and swallow ..." Need to combine "and swallow" with previous "chew"
                            ((AND_CONJ|and) \s)?
                            (?P<second_directive>DIRECTIVE)
                            (?![a-z])
                        ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        left_context_pattern = re.compile('DIRECTIVE')
        prior_directive_obj = left_context_pattern.search(left_context)
        if not prior_directive_obj:
            return None
        prior_directive_start = prior_directive_obj.start()
        prior_directive = parse.position2struc(prior_directive_start)

        groupname2group = trim_dictionary(match_object.groupdict())
        second_directive_start = match_object.start('second_directive')
        second_directive = parse.position2struc(second_directive_start)

        if prior_directive.value == 'chew' and second_directive.value == 'swallow':
            prior_directive.value = 'chew and swallow'
            prior_directive.rules_used.append('broken_up_directives')
            return []
        elif prior_directive.value in ('mix', 'dissolve') and second_directive.value in ('take', 'give'):
            # E.g. "dissolve one capful in water and take by mouth one time daily"
            return []
        elif prior_directive.value == 'use' and second_directive.value in ('test'):
            prior_directive.value = second_directive.value
            prior_directive.rules_used.append('broken_up_directives')
            return []
        elif prior_directive.value in ('use', 'give', 'take') and second_directive.value in ('swish and swallow'):
            prior_directive.value = second_directive.value
            prior_directive.rules_used.append('broken_up_directives')
            return []
        elif prior_directive.value in ('use', 'give', 'take') and second_directive.value in ('swish and spit'):
            prior_directive.value = second_directive.value
            prior_directive.rules_used.append('broken_up_directives')
            return []
        elif prior_directive.value in ('use', 'give', 'take') and second_directive.value in ('swish'):
            prior_directive.value = second_directive.value
            prior_directive.rules_used.append('broken_up_directives')
            return []
        elif prior_directive.value in ('swish') and second_directive.value in ('swallow'):
            prior_directive.value = 'swish and swallow'
            prior_directive.rules_used.append('broken_up_directives')
            return []
        elif prior_directive.value in ('swish') and second_directive.value in ('spit'):
            prior_directive.value = 'swish and spit'
            prior_directive.rules_used.append('broken_up_directives')
            return []

        return None




    rule = Rule_ExtractStrucs(   name = 'broken_up_directives',
                                    search_patterns = [pattern0],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule

def rule_assembly_dir_attribute():

    pattern0 = re.compile(r'''
                            (?<![a-z])
                            (?P<attribute>
                                sparingly               |
                                (?P<thinly>
                                    thinly                                  |
                                    (in \s)? (a \s)? thin \s (layer|film)
                                )                       |
                                lightly                 |
                                liberally               |
                                (?P<together>together)  |
                                (?P<first_dose>
                                    (as|for|in) \s
                                    (a \s|your \s)?
                                    (the \s)?
                                    (?P<quant2>QUANT)    # Quant better be = 1  and of type ordinal
                                    \s?
                                    day(\'?s)?          # "first day's dose"
                                    \s dose
                                )                       |
                                (?P<as_single_dose>         # could be as_single_dose or first_dose, depending if Quant is ordinal or not
                                    (as|for|in) \s
                                    (a \s)?
                                    (
                                        single          |
                                        (?P<quant>QUANT)    # Quant better be = 1 (and not "first dose" if it to mean "single dose")
                                    )
                                    \s dose
                                )
                            )
                            (?![a-z])
                        ''', re.X)

    pattern1 = re.compile(r'''                      #  "all at once" or "all at one time", which means "as single dose"
                            (?<![a-z])
                            (?P<as_single_dose>
                                (all \s)?
                                at \s
                                (
                                    once                        |
                                    (?P<quant>QUANT) \s time        # E.g. all at one time (Quant better be = 1)
                                )
                            )
                            (?![a-z])
                        ''', re.X)

    pattern2 = re.compile(r'''                      #  this is your first dose or "for the first dose"
                            (?<![a-z])
                            (?P<first_dose>
                                (
                                    (as \s)? your   |
                                    for \s the
                                )
                                \s
                                (?P<quant2>QUANT)    # Quant better be = 1  and of type ordinal
                                \s
                                dose
                            )
                            (?![a-z])
                        ''', re.X)

    pattern3 = re.compile(r'''                      #  take this once
                            (?<![a-z])
                            (?P<once_only>
                                (?P<quant>QUANT)    # Quant better be = 1 and of type INT
                                \s
                                time
                                (\s only)?
                                \.?
                            )
                            (?![a-z])
                        ''', re.X)


    pattern4 = re.compile(r'''                      #  inject per sliding scale
                            (?<![a-z])
                            (?P<per_sliding_scale>
                                (
                                   (
                                        per?            |
                                        on (\s a)?       |
                                        for              |
                                        using (\s the)?
                                   )
                                   \s
                                )?
                                sliding \s scales?
                            )
                            (?![a-z])
                        ''', re.X)

    pattern5 = re.compile(r'''                      #  apply a small amount
                            (?<![a-z])
                            (?P<small_amount>
                               (a \s)?
                               small \s amount s?
                            )
                            (?![a-z])
                        ''', re.X)


    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())
        attribute = match_object.group()

        if 'per_sliding_scale' in groupname2group:  # pattern 4: inject per sliding scale
            attribute = 'per_sliding_scale'
        elif 'DIRECTIVE' not in left_context and 'ANAPHORA' not in left_context:  # If there is not directive to the left, better then be "this is your first dose" (from Dictionary)
            return None

        if 'thinly' in groupname2group:
            attribute = 'thinly'

        if 'together' in groupname2group:               # "together" means "as one dose" only if it closely follows DOSE
                                                        # e.g. "take 2 tabs by mouth together
            required_left_context = re.compile('DOSE\)?(\sROUTE)?\s$')
            if not required_left_context.search(left_context):
                return None
            else:
                attribute = 'together'

        if 'as_single_dose' in groupname2group:             # pattern0 or pattern1
            if 'quant' in groupname2group:                  # Quant better be = 1
                quant_start = match_object.start('quant')
                quant = parse.position2struc(quant_start)
                if quant.value != 1:
                    return None
                if quant.num_type == 'int':                 # 'as one dose', which means "single dose"
                    attribute = 'as_single_dose'
                elif quant.num_type == 'ord':               # 'as first dose'
                    attribute = 'first_dose'
                else:
                    return None
            else:
                attribute = 'as_single_dose'

        if 'first_dose' in groupname2group:             # pattern2: this is your first dose
            quant_start = match_object.start('quant2')  # Quant better be = 1
            quant = parse.position2struc(quant_start)
            if quant.value != 1 or quant.num_type != 'ord':
                return None
            attribute = 'first_dose'



        if 'once_only' in groupname2group:      # pattern3 -- take this medicine once.
            quant_start = match_object.start('quant')
            quant = parse.position2struc(quant_start)
            if quant.value != 1 or quant.num_type != 'int':
                return None
            good_right_context_pattern = re.compile('\s*$|^\s*(CALENDAR_EVENT|TIMING|DIR_ATTRIBUTE)') # Example of dir_attribute on the right: "take 1 tablet once for 1 dose"
            bad_left_context_pattern = re.compile('AND_CONJ\s?$')       # omit things such as "take 3 times during the day and once at bedtime" because "once" here doesn't mean "one time only"
            if bad_left_context_pattern.search(left_context) or not good_right_context_pattern.match(right_context):
                return None
            attribute = 'once_only'


        if 'small_amount' in groupname2group:           # pattern5 apply a small amount
            attribute = 'small_amount'

        dir_attribute = DirAttribute(value = attribute)

        return [dir_attribute]

    rule = Rule_ExtractStrucs(   name = 'dir_attribute',
                                    search_patterns = [pattern0, pattern1, pattern2, pattern3, pattern4, pattern5],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule


def rule_assembly_range_identifcation():

    pattern = re.compile(r'''
                            (?P<low>QUANT)
                            \s?
                            (?P<range_type>to|or|\-|thru|through)
                            (\s?(to|or|\-|thru|through) )?          # sometimes they use double ranges "one -to two"
                            \s?
                            (?P<high>QUANT)
                        ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        bad_left_context_pattern = re.compile('QUANT\s*(\-|\/)*\s*$')
        bad_right_context_pattern = re.compile('^\s*(\-|\/)*\s*QUANT')
        if bad_left_context_pattern.search(left_context) or bad_right_context_pattern.search(right_context):
            return None
        groupname2group = trim_dictionary(match_object.groupdict())
        range_type = groupname2group['range_type']
        low_start = match_object.start('low')
        low = parse.position2struc(low_start)
        high_start = match_object.start('high')
        high = parse.position2struc(high_start)
        if low.num_type == 'range' or high.num_type == 'range':
            return None
        if low.num_type != 'var' and high.num_type != 'var' and low.value >= high.value:
            if range_type == '-' and low.num_type == 'int' and high.num_type == 'frac' and high.integer == 0:   # Sometimes a fraction such as 1.5 is represented as "1-1/2", so then fix it.
                high.integer = low.value
                return [high]
            else:
                return None
        if range_type == 'or':
            range_type = 'or'
        else:
            range_type = 'to'
        range_struc = RangeQuant(range_type = range_type, low = low, high = high)
        range_struc.constituents = [low, range_type, high]

        return [range_struc]

    rule = Rule_ExtractStrucs(   name = 'numeric range',
                                    search_patterns = [pattern],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule

def rule_assembly_dose():

    pattern0 = re.compile(r'''
                            (               # optional: there is a bizarre but frequent occurence of " 0 5ml". This should be 0.5ml
                                (?<=\s)
                                (?P<preceding_quant>QUANT)
                                \s
                            )?
                            (?P<quant>QUANT)
                            (
                                \s of (\s an?)? \s      |   # e.g. "take half of a tab"
                                \s a \s                 |   # e.g. "take half a teaspoon"
                                \s? \-? \s?                 # e.g. "take 1-tablets"
                            )
                            (
                                (?P<form>FORM) (\s dose)?  |     # e.g. "use 1 15mg dose"
                                (?P<doses>doses?)
                            )
                        ''', re.X)

    pattern1 = re.compile(r'''          # up to 3 tablets.
                            (?P<up_to>(\-\s?)? up \s? to \s)
                            (?P<dose>DOSE)
                        ''', re.X)

    pattern2 = re.compile(r'''          # with word "contents", eg. "inhale contents of 1 capsule"
                            (?P<contents_of>(the \s)? contents? \s of \s)
                            (
                                (?P<dose_contents_of>DOSE)      |
                                (?P<form_contents_of>FORM)
                            )
                        ''', re.X)


    pattern3 = re.compile(r'''          # Dose range expressed in the form "dose to dose", eg "take 1 tablet to 2 tablets daily"
                            (?P<dose_low>DOSE)
                            \s?
                            (to|\-)
                            \s?
                            (?P<dose_high>DOSE)
                        ''', re.X)

    pattern4 = re.compile(r'''          # Special pattern for injections and administration of weighted doses, e.g. "inject one 10mcg dose"
                            (?P<quant_one>QUANT) \s
                            (?P<quant>QUANT)
                            (
                                \s of (\s an?)? \s      |   # e.g. "take half of a tab"
                                \s a \s                 |   # e.g. "take half a teaspoon"
                                \s? \-? \s?                 # e.g. "take 1-tablets"
                            )
                            (?P<form>FORM)
                            (\s dose)?
                        ''', re.X)

    pattern5 = re.compile(r'''          # Insert missing (implied) quant = 1 to create a dose
                                        # E.g.  "take tablet by mouth as directed by physician", "insert ring vaginally for 3 weeks"
                            (?P<implied_quant>
                                (?P<directive>DIRECTIVE)
                                \s
                                (
                                    (an|a|the)
                                    \s
                                )?
                                (?P<form>FORM)
                            )
                        ''', re.X)




    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        # Avoid classifying "taper 5-4-3-2-1 tabs daily" as QUANT-QUANT-QUANT-DOSE with dose.quant.value = 1
        bad_left_context_pattern = re.compile('QUANT\s*(\-|\/)*\s*$')
        if bad_left_context_pattern.search(left_context):
            return None
        groupname2group = trim_dictionary(match_object.groupdict())


        if 'quant' in groupname2group:      # pattern0
            quant_start = match_object.start('quant')
            quant = parse.position2struc(quant_start)
            if 'preceding_quant' in groupname2group:
                preceding_quant_start = match_object.start('preceding_quant')
                preceding_quant = parse.position2struc(preceding_quant_start)
                if preceding_quant.num_type == 'int' and preceding_quant.value == 0 and quant.num_type == 'int' and quant.value < 10 and quant.value > 0:
                    quant.num_type = 'decimal'
                    quant.value = quant.value / 10.0
                    quant.value_repr = str(quant.value)


        if 'form' in groupname2group:
            form_start = match_object.start('form')
            form = parse.position2struc(form_start)

        if 'doses' in groupname2group:
            return None
            #if quant.num_type == 'ord':     # this is "first_dose", which is not a form but a directive attribute
            #    return None
            #form_value = groupname2group['doses']
            #if form_value[-1] == 's':
            #    plurality = 'plural'
            #else:
            #    plurality = 'singular'
            #form = Form(form_name = 'dose', plurality = plurality, constituents = [form_value])
            #form.rules_used.append('dose')
        elif 'up_to' in groupname2group:        # pattern1
            dose_start = match_object.start('dose')
            dose = parse.position2struc(dose_start)
            dose_value = dose.quant.value
            if dose.quant.num_type not in ('int', 'var') or dose_value <= 1:
                return None
            dose.up_to = True
            return [dose]
        elif 'contents_of' in groupname2group:        # pattern2  Could be "inhale contents of 1 capsule" or "inhale content capsule"
            if 'dose_contents_of' in groupname2group:
                dose_start = match_object.start('dose_contents_of')
                dose = parse.position2struc(dose_start)
            else:       # it's FORM, e.g. "inhale contents capsule"
                form_start = match_object.start('form_contents_of')
                form = parse.position2struc(form_start)
                quant = Quant(num_type = 'int', value = 1)
                dose = Dose(quant = quant, form = form)
            return [dose]
        elif 'dose_low' in groupname2group:        # pattern3: Dose to Dose
            dose_low_start = match_object.start('dose_low')
            dose_low = parse.position2struc(dose_low_start)
            dose_high_start = match_object.start('dose_high')
            dose_high = parse.position2struc(dose_high_start)
            if dose_low.form.value == dose_high.form.value and dose_low.quant.value < dose_high.quant.value:
                quant = RangeQuant(range_type = 'to', low = dose_low.quant, high = dose_high.quant)
                form = dose_low.form
            else:
                return None
        elif 'quant_one' in groupname2group:    # pattern such as "inject one 10mcg dose". the quant_one better be 1 and the form better be weight/volume based
            quant_one_start = match_object.start('quant_one')
            quant_one = parse.position2struc(quant_one_start)
            if quant_one.value != 1 or form.value not in ('cc', 'gram', 'mcg', 'mg', 'ml', 'ounce'):
                return None

        elif 'implied_quant' in groupname2group:    # pattern5 such as "take tablet by mouth as directed by physician", "insert ring vaginally for 3 weeks"
            directive_start = match_object.start('directive')
            directive = parse.position2struc(directive_start)
            if form.value in ('cc', 'gram', 'mcg', 'mg', 'ml', 'ounce'):
                return None
            elif form.plurality != 'singular':
                return None
            quant = Quant(num_type = 'int', value = 1)
            space = Struc(' ', [])
            dose = Dose(quant = quant, form = form)
            return [directive, space, dose]



        dose = Dose(quant = quant, form = form)
        return [dose]

    rule = Rule_ExtractStrucs(   name = 'dose',
                                    search_patterns = [pattern0, pattern1, pattern2, pattern3, pattern4, pattern5],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule



def rule_assembly_dose_omitted_form():
    """ Process doses with implied form. E.g. "Take 2 tabs in the AM and 3 in the PM". Here "3" means "3 tabs"
    """

    pattern = re.compile(r'''
                            (?P<quant>QUANT)
                            (?=
                                \s?
                                (
                                    ROUTE           |
                                    FREQ            |
                                    PERIODICITY     |
                                    TIMING          |
                                    CALENDAR_EVENT
                                )
                            )
                        ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        def get_first_known_dose_form_prior_to_this_struc(struc_position):
            """ Returns the Form() struc that occurs in prior to this struc_position in the parse.strucs list, closest to the position.
                If none are found, returns the best compatible with verb form.
            """

            doses = [struc for struc in parse.strucs[:struc_position] if struc.label == 'DOSE']
            if doses:
                last_dose = doses[-1]
                form = last_dose.form
                return form

            forms = [struc for struc in parse.strucs[:struc_position] if struc.label == 'FORM']
            if forms:
                return forms[-1]


            directives  = [struc for struc in parse.strucs[:struc_position] if struc.label == 'DIRECTIVE']
            for directive in directives:
                form_name = Directive.directive_to_likely_form.get(directive.value)
                if form_name:
                    form = Form(form_name = form_name, plurality = '', constituents = [])
                    return form

            return None
        groupname2group = trim_dictionary(match_object.groupdict())
        quant_start = match_object.start('quant')
        quant = parse.position2struc(quant_start)

        quant_index = parse.strucs.index(quant)
        previous_form = get_first_known_dose_form_prior_to_this_struc(quant_index)
        if previous_form:
            form = Form(form_name = previous_form.value, plurality = '', constituents = [])
            form.rules_used.append('*form_deduced_from_quant_in_struc_id_dose_omitted_form_deduced()')
            if quant.value != 1:
                form.plurality = 'plural'
            else:
                form.plurality = 'singular'
            dose = Dose(quant, form)
            return [dose]
        else:
            return None

    rule = Rule_ExtractStrucs(   name = 'dose_omitted_form_deduced',
                                    search_patterns = [pattern],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule


def rule_remove_parenthesis_duplicate_dose():
    """ Deal with cases such as "take 1 teaspoonful (5ml) by mouth on day 1...", "take 2 tablets (1 gram) together ", "take 1/2 tablet (=25mg dose) by mouth every day.", etc.
    """

    pattern0 = re.compile(r'''              # e.g "take one and 1/2 tablets (1.5 tablets)"
                            (?<!MAX)        # don't do it for MAXDOE (DOSE)
                            (?P<dose1> DOSE)
                            \s?
                            \(
                            \=? \s?
                            (?P<dose2> DOSE)
                            (\s dose)?
                            \s?
                            \)
                        ''', re.X)

    pattern1 = re.compile(r'''              # remove dictionary items, eg. "mix 2 teaspoons (1 teaspoon = 5 ml)"
                            (?<!MAX)        # don't do it for MAXDOE (DOSE)
                            (?P<dose1> DOSE)
                            \s?
                            \(
                            (?P<dose2> DOSE) \s?
                            \=? \s?
                            (?P<dose3> DOSE)
                            \s?
                            \)
                        ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        dose1_start = match_object.start('dose1')
        dose1 = parse.position2struc(dose1_start)

        return [dose1]


    rule = Rule_ExtractStrucs(  name = 'remove_parenthesis_duplicate_dose',
                                search_patterns = [pattern0, pattern1],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)

    return rule

def rule_assembly_substrate():

    pattern = re.compile(r'''
                            (?<![a-z])
                            (with | in(to)? ) \s        # mandatory
                            (
                                (?P<dose> DOSE)         # Dose is optional
                                (\s of)?
                                \s?
                            )?
                            ( (a \s)? (full \s)? glass \s (of \s)? )?       # e.g. Mix 1 capful in 8 oz glass of water" or "dissolve 1 packet in a full glass of water"
                            (?P<value>                  # mandatory
                                (?P<water>
                                    water \s or \s juice    |
                                    juice \s or \s water    |
                                    (plain \s)? water
                                )                                   |
                                (?P<gatorade>
                                    (
                                        (apple \s)?
                                        juice
                                        \s
                                        (or \s)?
                                    )?
                                    (
                                        gai?t.rade?         |
                                        gatorate            |
                                        gatarate            |
                                        powerade
                                    )
                                )                                   |
                                (?P<liquid>
                                    (any \s)? (clear \s)?
                                    (
                                        liquid              |
                                        fluid               |
                                        beverage            |
                                        juice
                                    )
                                )                                   |
                                (?P<diluent>diluent)
                            )
                            (?![a-z])
                        ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())
        if 'dose' in groupname2group:
            dose_start = match_object.start('dose')
            dose = parse.position2struc(dose_start)
        else:
            dose = None

        if 'water' in groupname2group:
            value = 'water'
        elif 'gatorade' in groupname2group:
            value = 'gatorade'
        elif 'diluent' in groupname2group:
            value = 'diluent'
        else:
            value = 'liquid'

        substrate = Substrate(value, dose)

        return [substrate]

    rule = Rule_ExtractStrucs(   name = 'substrate',
                                    search_patterns = [pattern],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule



def rule_assembly_spray_n_times():
    """ "Spray N times" is equivalent to "Use N sprays". I.e. we convert "n times" (which looks like freq) to a dose.
    """

    pattern = re.compile(r'''
                            (?P<directive>DIRECTIVE)
                            \s
                            (?P<quant>QUANT)
                            (\s|\-)?
                            times?
                            (?![a-z])
                            ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        groupname2group = trim_dictionary(match_object.groupdict())

        quant_start = match_object.start('quant')
        quant = parse.position2struc(quant_start)

        directive_start = match_object.start('directive')
        directive = parse.position2struc(directive_start)

        if directive.value != 'spray':
            return None
        if quant.value == 'range' or quant.value > 1:
            plurality = 'plural'
        else:
            plurality = 'singular'
        form = Form(form_name = 'spray', plurality = plurality, constituents = [quant.value_repr, 'times'])
        dose = Dose(quant = quant, form = form)
        directive = Directive(value = 'use')
        directive.constituents = [directive.value]

        space = Struc(label = ' ')
        return [directive, space, dose]


    rule = Rule_ExtractStrucs(   name = 'spray n times -> use n sprays',
                                    search_patterns = [pattern],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule









def rule_assembly_time_unit():

    pattern0 = re.compile(r'''
                            (?<![a-z])
                            (?P<time_unit>
                                hour|hr|day|week|wk|month|minute|second
                            )
                            (?P<plural>
                                (s|\(s\))
                            )?
                            (?![a-z])
                            ''', re.X)

    pattern1 = re.compile(r'''                 # Special case for minutes. If abbreviated "min" we disambiguate it as minutes if preceded by a QUANT
                            (?<=QUANT)
                            (\s|\-)?
                            (?P<time_unit>
                                min
                            )
                            (?P<plural>
                                (s|\(s\))
                            )?
                            (?![a-z])
                            ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        groupname2group = trim_dictionary(match_object.groupdict())
        time_unit = match_object.group('time_unit')
        if time_unit in ('hour', 'hr'):
            time_unit = 'hour'
        elif time_unit in ('week', 'wk'):
            time_unit = 'week'
        elif time_unit in ('min', 'minute'):
            time_unit = 'minute'
        else:
            time_unit = time_unit.lower()

        if 'plural' in groupname2group:
            plurality = groupname2group['plural']
            if plurality == 's':
                plurality = 'plural'
            else:
                plurality = 'plurality_either'
        else:
            plurality = 'singular'

        time_unit = TimeUnit(value = time_unit, constituents = [match_object.group()], plurality = plurality)

        return [time_unit]

    rule = Rule_ExtractStrucs(   name = 'time unit',
                                    search_patterns = [pattern0, pattern1],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule


def rule_assembly_freq():

    pattern0 = re.compile(r'''
                            (   # covers strings such as "2 times/wk", "twice daily", "twice every day"
                                (?P<quant>QUANT)
                                \s? (\-\s?)?                         # could be "2-times a day"
                                times? \s?
                                (
                                    an?     |
                                    per     |
                                    \/
                                )?
                                \s?
                                (   # Option 1: "2 times an hour" or "twice every day"
                                        (
                                            an?                 |
                                            per                 |
                                            \/                  |
                                            (?P<every>every|each)
                                        )
                                        \s?
                                        (?P<time_unit1>TIMEUNIT)
                                |   # Option 2: "daily"
                                        (?P<time_unit2>hourly|daily|weekly|monthly)   # covers strings such as "5 times weekly"
                                )
                             |      # Option 3: covers implied freq of 1, e.g "take 1 tab by mouth daily"
                                (?P<time_unit3>hourly|daily|weekly|monthly)
                            )
                           ''', re.X)

    pattern1 = re.compile(r'''  # "once or twice a day" or "2 times to 3 times a day" where "3 times daily" was already processed as FREQ
                            (?P<quant_low>QUANT)
                            (\s|\s?\-\s?)?
                            times? \s?
                            (?P<range_quant>to|or|\-)         # mandatory
                            \s?
                            (?P<high_freq>FREQ)
                            ''', re.X)



    pattern2 = re.compile(r'''  # "2 times a day to 3 times a day" where "3 times daily" was already processed as FREQ
                            (?P<low_freq>FREQ)
                            (\s (?P<periodicity>PERIODICITY))?  # e.g. "2 times DAILY to 3 times daily"
                            \s?
                            (?P<range_freq>to|or|\-)         # mandatory
                            \s?
                            (?P<high_freq>FREQ)
                            ''', re.X)

    pattern3 = re.compile(r'''          # E.g., "up to 3 times daily".
                            (?P<up_to>(\-\s)? up \s? to \s)
                            (?P<freq>FREQ)
                        ''', re.X)

    pattern4 = re.compile(r'''          # E.g., "daily to 3 times daily". We change it to "1 to 3 times daily"
                            (?<=PERIODICITY\s)
                            (?P<implied_range>(\-\s)? to \s)
                            (?P<freq>FREQ)
                        ''', re.X)


    pattern5 = re.compile(r'''          # Bad use of "day" without indefinite article instead of "a day" or "daily" in the clear FREQ context. E.g., "use 3 times day.$"
                            (?P<quant>QUANT) \s
                            times \s
                            (?P<time_unit1>TIMEUNIT)
                            \s?\.? $
                        ''', re.X)


    pattern6 = re.compile(r'''          # Covers patterns with dayparts. E.g. "2 times a night", "once every night", "once every morning"
                            (?P<quant>QUANT) \s
                            times? \s
                            (
                                an?                 |
                                per                 |
                                \/                  |
                                (?P<every>every|each)
                            )
                            \s?
                            (?P<day_part>morning|evening|night|afternoon|bedtime)s?
                        ''', re.X)

    pattern7 = re.compile(r'''  # covers strings such as "2 days per wk" which means "2 times per week"
                                (?P<timeinterval>TIMEINTERVAL)
                                \s?
                                (
                                    an?     |
                                    per     |
                                    \/
                                )
                                \s?
                                (?P<time_unit1>TIMEUNIT)
                           ''', re.X)


    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())

        bad_left_context_pattern = re.compile('[a-zA-Z]$')
        bad_right_context_pattern = re.compile('^[a-zA-Z]')
        if bad_left_context_pattern.search(left_context) or bad_right_context_pattern.search(right_context):
            return None

        # There are cases when "daily" or "weekly" without a preceding quantity do NOT imply quantity 1.
        # This typically happens in the context of "max daily dose" type of expressions,
        # e.g. "max daily dose 3/24" where "daily" does not mean "1 every day"
        bad_left_context_pattern_for_implied_time_freq_1 = re.compile('max|exceed|mdd|not? more', re.I)

        if 'time_unit3' in groupname2group and bad_left_context_pattern_for_implied_time_freq_1.search(left_context):
            return None

        if 'quant' in groupname2group:
            quant_start = match_object.start('quant')
            quant = parse.position2struc(quant_start)
        else:
            quant = Quant(num_type = 'int', value = 1)

        if 'time_unit1' in groupname2group:
            time_unit_start = match_object.start('time_unit1')
            time_unit = parse.position2struc(time_unit_start)
        elif 'time_unit2' in groupname2group or 'time_unit3' in groupname2group:
            if 'time_unit2' in groupname2group:
                time_unit_txt = groupname2group['time_unit2']
            else:
                time_unit_txt = groupname2group['time_unit3']

            if time_unit_txt == 'hourly':
                time_unit = TimeUnit(value = 'hour', constituents = [time_unit_txt])
            elif time_unit_txt == 'daily':
                time_unit = TimeUnit(value = 'day', constituents = [time_unit_txt])
            elif time_unit_txt == 'weekly':
                time_unit = TimeUnit(value = 'week', constituents = [time_unit_txt])
            else:
                time_unit = TimeUnit(value = 'month', constituents = [time_unit_txt])
        elif 'day_part' in groupname2group:         # pattern6
            time_unit = TimeUnit(value = 'day', constituents = [], plurality = 'singular')
            freq = Freq(quant = quant, time_unit = time_unit)
            quant1 = Quant(num_type = 'int', value = 1)
            time_unit1 = TimeUnit(value = 'day', constituents = [], plurality = 'singular')
            periodicity = Periodicity(quant = quant1, time_unit = time_unit1)
            landmark = groupname2group['day_part']
            relation = 'at' if landmark == 'bedtime' else ''
            timing = Timing(landmark = landmark, relation = relation)
            space1 = Struc(label = ' ')
            space2 = Struc(label = ' ')
            return [freq, space1, periodicity, space2, timing]

        if 'range_quant' in groupname2group:        # pattern1
            # We are dealing with pattern1: eg. "2 times or 3 times daily" where "3 times daily" was already processed as FREQ
            # We need to create a new Quant object which is going to be the range
            range_type = groupname2group['range_quant']
            if range_type != 'or':
                range_type = 'to'
            quant_low_start = match_object.start('quant_low')
            quant_low = parse.position2struc(quant_low_start)
            high_freq_start = match_object.start('high_freq')
            high_freq = parse.position2struc(high_freq_start)
            quant_high = high_freq.quant
            time_unit = high_freq.time_unit
            quant = RangeQuant(range_type = range_type, low = quant_low, high = quant_high)
            quant.constituents = [quant_low, range_type, quant_high]

        if 'range_freq' in groupname2group:         # pattern2
            # We are dealing with pattern2: FREQ to FREQ. We need to verify that time_units are the same
            # Then we modify the freq.quant of the last FREQ
            range_type = groupname2group['range_freq']
            if range_type != 'or':
                range_type = 'to'
            low_freq_start = match_object.start('low_freq')
            low_freq = parse.position2struc(low_freq_start)
            high_freq_start = match_object.start('high_freq')
            high_freq = parse.position2struc(high_freq_start)
            if not low_freq.time_unit.value == high_freq.time_unit.value:
                return None
            time_unit = high_freq.time_unit
            quant_low_freq = low_freq.quant
            quant_high_freq = high_freq.quant
            if quant_low_freq.num_type != 'int' or quant_high_freq.num_type != 'int':
                return None
            low_value = quant_low_freq.value
            high_value = quant_high_freq.value
            if not low_value < high_value:
                return None
            if 'periodicity' in groupname2group: # check that Periodicity = "daily"
                periodicity_start = match_object.start('periodicity')
                periodicity = parse.position2struc(periodicity_start)
                if periodicity.quant.value != 1:
                    return None
            quant = RangeQuant(range_type = range_type, low = quant_low_freq, high = quant_high_freq)
            quant.constituents = [low_freq, range_type, high_freq]

        if 'timeinterval' in groupname2group:         # pattern7
            # 2 days per week
            timeinterval_start = match_object.start('timeinterval')
            timeinterval = parse.position2struc(timeinterval_start)
            if time_unit.value != 'week' or timeinterval.time_unit.value != 'day' or timeinterval.quant.num_type != 'int':
                return None
            quant = timeinterval.quant

        if 'up_to' in groupname2group:        # pattern3
            freq_start = match_object.start('freq')
            freq = parse.position2struc(freq_start)
            if freq.quant.num_type not in ('int', 'var') or freq.quant.value <= 1:
                return None
            freq.up_to = True
            return [freq]

        if 'implied_range' in groupname2group:        # pattern4
            freq_start = match_object.start('freq')
            freq = parse.position2struc(freq_start)
            time_unit = freq.time_unit
            freq_value = freq.quant.value
            if freq.quant.num_type != 'int' or freq_value <= 1:
                return None
            freq_low_range = Quant(num_type = 'int', value = 1)
            quant = RangeQuant(range_type = 'to', low = freq_low_range, high = freq.quant)

        freq = Freq(quant = quant, time_unit = time_unit)

        if 'every' in groupname2group or 'time_unit2' in groupname2group or 'time_unit3' in groupname2group:
            # There is both freq and implied periodicity, because we are dealing with cases such as "twice daily" or "twice every day"
            # which mean "take 2 tabs per day. Take this every day".
            quant_periodicity = Quant(num_type = 'int', value = 1)
            time_unit_periodicity = TimeUnit(value = time_unit.value, constituents = [])
            periodicity = Periodicity(quant = quant_periodicity, time_unit = time_unit_periodicity)
            space = Struc(label = ' ')
            if 'time_unit3' in groupname2group:
                # Although I waffle on it, currently I think that "take 2 tabs daily" does not imply FREQ of 1/day, only periodicity of "every 1 day".
                return [periodicity]
            return [freq, space, periodicity]
        else:
            return [freq]

    rule = Rule_ExtractStrucs(   name = 'freq',
                                    search_patterns = [pattern0, pattern1, pattern2, pattern3, pattern4, pattern5, pattern6, pattern7],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule


def rule_assembly_periodicity():

    pattern0 = re.compile(r'''
                                (?<![a-z])
                                ((?P<quant_once>QUANT) \s time \s)?     # "once every 6 hours", which is preprocessed as "1 time every 6 hours"
                                every \s?
                                (?P<time_interval>TIMEINTERVAL)
                                (?![a-z])
                            ''', re.X)

    pattern1 = re.compile(r'''          # covers cases such as "every other day", "every 3rd day", and "every day", "every morning", "every other night"
                                (?<![a-z])
                                (at \s)?                    # e.g. "at every bedtime"
                                (every|each) \s?             # mandatory
                                (                           # optional "other" or "2nd"
                                    (?P<other>other \s)         |
                                    (?P<ordinal>QUANT \s?)
                                )?
                                (
                                    (?P<timeunit>TIMEUNIT)                                  |
                                    (?P<day_part>morning|evening|night|afternoon|bedtime)   |
                                    (?P<timing>TIMING)
                                )s?
                                (?![a-z])
                        ''', re.X)

    pattern2 = re.compile(r'''          # covers cases such as "every Monday" or "each Sunday"
                                (?<![a-z])
                                (every|each)
                                (?=
                                    \s
                                    (?P<specific_day>SPECIFIC_DAY)
                                )
                        ''', re.X)


    pattern3 = re.compile(r'''          # covers cases without Quant, (implied 1) such as "take 2 tabs a day", "2 a day", "2 per week", take 2/day (which mean "take 2 tabs every day"
                                        # But can't be preceded by "times", which would be Freq
                                (?<![a-z])
                                (?P<implied_quant>  # mandatory
                                    \s a \s     |
                                    \s per \s   |
                                    \/ \s?
                                )
                                (?P<timeunit>TIMEUNIT)
                                (?![a-z])
                        ''', re.X)


    pattern4 = re.compile(r'''          # E.g., "up to every 4 hours".
                            (?P<up_to>(\-\s)? up \s? to \s)
                            (?P<max_period>PERIODICITY)
                        ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())

        if 'time_interval' in groupname2group:      # pattern0
            time_interval_start = match_object.start('time_interval')
            time_interval = parse.position2struc(time_interval_start)
            if 'quant_once' in groupname2group:
                quant_once_start = match_object.start('quant_once')
                quant_once = parse.position2struc(quant_once_start)
                if not quant_once.value == 1:
                    return None
            quant = time_interval.quant
            time_unit = time_interval.time_unit
            if quant.num_type == 'int':
                plurality = 'plural' if quant.value > 1 else 'singular'
                time_unit.plurality = plurality
        elif 'timeunit' in groupname2group or 'day_part' in groupname2group or 'timing' in groupname2group:       # pattern1: no explicit TIMEINTERVAL is given
            if 'other' in groupname2group:      # "every other day" means "every second day"
                quant = Quant(num_type = 'int', value = 2)
            elif 'ordinal' in groupname2group:
                ordinal_start = match_object.start('ordinal')
                ordinal = parse.position2struc(ordinal_start)
                quant = ordinal
                quant.num_type == 'int'  # change num_type to integer, because "every third week" = "every 3 weeks"
            else:                               # "every week" means "every 1 week"
                quant = Quant(num_type = 'int', value = 1)

            plurality = 'plural' if quant.value > 1 else 'singular'

            if 'timeunit' in groupname2group:
                time_unit_start = match_object.start('timeunit')
                time_unit = parse.position2struc(time_unit_start)
                if quant.num_type == 'int':
                    time_unit.plurality = plurality # modify time_unit plurality to fit the actual number seen

            else:       # 'day_part', e.g. "every morning" or "every other night" or Timing (i.e. it is already processed)
                if 'day_part' in groupname2group:
                    landmark = groupname2group['day_part']
                    relation = 'at' if landmark == 'bedtime' else ''
                    timing = Timing(landmark = landmark, relation = relation)
                elif 'timing' in groupname2group:
                    timing_start = match_object.start('timing')
                    timing = parse.position2struc(timing_start)
                    landmark = timing.landmark
                    if landmark not in ('morning', 'evening', 'night', 'afternoon', 'bedtime', 'lunch'):
                        return None

                time_unit = TimeUnit(value = 'day', constituents = [], plurality = plurality)
                periodicity = Periodicity(quant = quant, time_unit = time_unit)

                space = Struc(label = ' ')
                return [periodicity, space, timing]


        elif 'specific_day' in groupname2group:       # pattern2: e.g. "every Monday"
            specific_day_start = match_object.start('specific_day')
            specific_day = parse.position2struc(specific_day_start)
            if specific_day.typ == 'day_of_week':
                time_unit = TimeUnit(value = 'week', constituents = [], plurality = 'singular')
                quant = Quant(num_type = 'int', value = 1)
            else:
                return None

        elif 'up_to' in groupname2group:        # pattern4
            periodicity_start = match_object.start('max_period')
            periodicity = parse.position2struc(periodicity_start)
            if periodicity.quant.num_type not in ('int', 'var') or periodicity.quant.value <= 1:
                return None
            periodicity.up_to = True
            return [periodicity]

        if 'implied_quant' in groupname2group:      # pattern3: "2 tabs a day"
            bad_left_context1 = re.search("(times?|FREQ|TIMEINTERVAL)\s?$", left_context)    # could be part of FREQ (2 days per week) or of MaxDose
            bad_left_context2 = re.search("max|mdd|more|exceed", left_context)  # could be MaxDose
            if bad_left_context1 or bad_left_context2:
                return None
            quant = Quant(num_type = 'int', value = 1)
            space = Struc(label = ' ')
            periodicity = Periodicity(quant = quant, time_unit = time_unit)
            return [space, periodicity]



        periodicity = Periodicity(quant = quant, time_unit = time_unit)

        return [periodicity]

    rule = Rule_ExtractStrucs(   name = 'periodicity',
                                    search_patterns = [pattern0, pattern1, pattern2, pattern3, pattern4],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule



def rule_assembly_route():

    pattern0 = re.compile(r'''
                            (?<![a-z])
                            (ROUTE \s)?                 # there are shocking number of cases of repeated routes. We just omit the first instance.
                            (?P<route>
                                (?P<externally>
                                    externally              |
                                    externaly               |
                                    \(? external \s (use    |
                                    (?P<external_use_dir>DIRECTIVE)) \s only \)?
                                )                                                   |
                                (?P<intramuscularly>
                                    intramuscularly                 |
                                    intramuscular (\s injection)?   |
                                    in(to)? (\s the)? \s muscles?
                                )                                                   |
                                (?P<intranasally>
                                    intranasall?y                   |
                                    intranasal                      |
                                    ((by|via) \s (the \s)?)? nasal \s route |
                                    in \s nostril                   |
                                    nasally                         |
                                    in \s (the \s)? nose
                                )                                                   |
                                (?P<orally>
                                    by \s mouth             |
                                    in \s (the \s)? mouth   |       # but this could also be the site "rinse in the mouth", even though it is often "dissolve
                                    oral(ly)?
                                )                                                   |
                                (?P<rectally>
                                    rectally            |
                                    (?<!to\s) rectal    |       # avoid "apply to rectal area", which is NOT rectal route but topical
                                    (into|in|to|via|per) \s (the \s)? rectum
                                )                                                   |
                                (?P<subcutaneously>
                                    subcutaneously          |
                                    subscutaneously         |
                                    subcutaneous            |
                                    sub\-?q                 |
                                    under (\s the)? \s skin |
                                    below (\s the)? \s skin
                                )                                                   |
                                (?P<topically>
                                    topical(ly)?            |
                                    locally
                                )                                                   |
                                (?P<sublingually>
                                    sublinguall?y           |
                                    sublingual              |
                                    under \s (the \s)? tongue
                                )                                                   |
                                (?P<transdermally>
                                    transdermall?y           |
                                    transdermal (\s route)?
                                )                                                   |
                                (?P<vaginally>
                                    (per|inside|into|intra|to|in) \s (the \s)? vaginally    |
                                    (inside|into|in) \s (the \s)? vagina                    |
                                    (in \s)? (by|per) \s (the \s)? vaginal \s route?        |
                                    (intra)? vaginall?y                                     |
                                    (?<!for\s) vaginal (?!\sopening)(?!\sitch)
                                )
                            )
                            (?![a-z])
                         ''', re.X)

    pattern1 = re.compile(r'''
                            (?<![a-z])
                            (by \s)?                    # e.g "by oral route"
                            (?P<prev_route>ROUTE)
                            (\s route)
                            (?![a-z])
                         ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())

        if 'externally' in groupname2group:
            route = 'externally'
        elif 'intranasally' in groupname2group:
            route = 'intranasally'
        elif 'intramuscularly' in groupname2group:
            route = 'intramuscularly'
        elif 'orally' in groupname2group:
            # make sure that the phrase (e.g. "in mouth" or "in the mouth") it doesn't refer to
            # the site = mouth, e.g. "rinse in mouth", which is not at all oral route.
            # There are even frequent uses of "rinse by mouth", which still means site = mouth.
            if 'rinse' in left_context or 'swish' in left_context:
                return None
            else:
                left_context_pattern = re.compile('(?P<left_dir>DIRECTIVE)')
                left_dir_match_obj = left_context_pattern.search(left_context)
                if left_dir_match_obj:
                    left_dir_start = left_dir_match_obj.start('left_dir')
                    left_dir = parse.position2struc(left_dir_start)
                    if left_dir.value in ('rinse', 'swish'):
                        return None
            route = 'orally'
        elif 'rectally' in groupname2group:
            route = 'rectally'
        elif 'subcutaneously' in groupname2group:
            route = 'subcutaneously'
        elif 'sublingually' in groupname2group:
            route = 'sublingually'
        elif 'topically' in groupname2group:
            route = 'topically'
        elif 'transdermally' in groupname2group:
            route = 'transdermally'
        elif 'vaginally' in groupname2group:
            route = 'vaginally'
        elif 'prev_route' in groupname2group:     # pattern1
            prev_route_start = match_object.start('prev_route')
            prev_route = parse.position2struc(prev_route_start)
            return [prev_route]



        route = Route(value = route, constituents = [match_object.group()])

        return [route]

    rule = Rule_ExtractStrucs(   name = 'route',
                                    search_patterns = [pattern0, pattern1],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule


def rule_assembly_timeinterval():

    pattern = re.compile(r'''
                        (?P<quant>QUANT)
                        \s? \-? \s?
                        (consecutive \s)?
                        (?P<time_unit>TIMEUNIT)
                        ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        quant_start = match_object.start('quant')
        quant = parse.position2struc(quant_start)
        if quant.num_type == 'ord':     # avoid ordinals, e.g. "first day" -- they are not intervals but Calendar_Events
            return None
        time_unit_start = match_object.start('time_unit')
        time_unit = parse.position2struc(time_unit_start)
        time_interval = TimeInterval(quant, time_unit)
        return [time_interval]

    rule = Rule_ExtractStrucs(   name = 'time interval',
                                    search_patterns = [pattern],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule


def rule_assembly_duration():
    """
    Duration period: how many days/hours/weeks to take the medication. E.g. "for 3 days", "for 4 more days", "for the next 3 days"
    """
    pattern0 = re.compile(r'''                  # for (the next)? 4 days
                            (?<![a-z])
                            for \s?             # mandatory
                            (?P<start_next>(the \s)? next \s)?
                            ((the \s)? (?P<start_first>QUANT) \s)?
                            (?P<additional>additional \s)?
                            (?P<time_interval>TIMEINTERVAL)
                            (?P<thereafter>\s
                                (
                                    there \s? after     |
                                    after \s that       |
                                    more
                                )
                            )?
                            (?![a-z])
                        ''', re.X)

    pattern1 = re.compile(r'''                  # for 4 more days
                            (?<![a-z])
                            (there \s? after \s)?       # e.g. "thereafter for 4 more days"
                            for \s
                            (?P<quant>QUANT) \s
                            (?P<additional>more|additional) \s         # mandatory
                            (?P<time_unit>TIMEUNIT)
                            (?![a-z])
                        ''', re.X)


    pattern2 = re.compile(r'''                  # for 4 days more
                            (?<![a-z])
                            (there \s? after \s)?       # e.g. "thereafter for 4 days more"
                            for \s
                            (?P<quant>QUANT) \s
                            (?P<time_unit>TIMEUNIT)
                            \s (?P<additional>more|additional)    # mandatory
                            (?![a-z])
                        ''', re.X)

    pattern3 = re.compile(r'''                  # for 4 NIGHTS (which means days)
                            (?<![a-z])
                            for \s
                            (?P<quant>QUANT) \s
                            (?P<start_next>(the \s)? next \s)?
                            (?P<additional>additional \s)?
                            (?P<nights>night(?P<plural>s)?)
                            (?![a-z])
                        ''', re.X)

    pattern4 = re.compile(r'''                  # on each of the next 4 days
                            (?<![a-z])
                            on \s each \s of \s                 # mandatory
                            (?P<start_next>(the \s)? next \s)   # mandatory
                            (?P<time_interval>TIMEINTERVAL)
                            (?P<thereafter>\s
                                (
                                    there \s? after     |
                                    after \s that
                                )
                            )?
                            (?![a-z])
                        ''', re.X)

    pattern5 = re.compile(r'''                  # bizarre but frequent occurence: at the end just says "3 days", e.g. 'take 1 tab daily 3 days"
                            (?P<terminal_time_interval>
                                (?P<time_interval>TIMEINTERVAL)
                                \s?\.?$
                            )
                        ''', re.X)

    pattern6 = re.compile(r'''                  # "up to 12 hours a day", "apply to affected area for up to 12 hours"
                            (?<![a-z])
                            (for \s)?
                            (?P<up_to>   # mandatory
                                up \s? to \s
                            )
                            (?P<time_interval>TIMEINTERVAL)
                            (?!\s(after|before))    # don't trap timing offset, e.g. "take up to 3 hours before lunch"
                            (?![a-z])
                        ''', re.X)


    pattern7 = re.compile(r'''                  # "Remove after 3 weeks" means "do whatever primary directive says for 3 weeks, then Remove", i.e. turn this to Duration Directive/remove
                                                # Same with "stop in 3 weeks"
                            (?<![a-z])
                            (?P<directive_remove>DIRECTIVE)
                            \s
                            (in | after)
                            \s
                            (?P<time_interval>TIMEINTERVAL)
                            (?![a-z])
                        ''', re.X)

    pattern8 = re.compile(r'''                  # "skip 7 days" or "skip a week" = "stop for 7 days", etc.
                            (?<![a-z])
                            (?P<skip>skip)
                            (\s (a|for|the \s next|the|next) )?
                            \s
                            (
                                (?P<time_interval>TIMEINTERVAL)         |
                                (?P<time_unit_no_quant>TIMEUNIT)        # means timeinterval with implied quant = 1
                            )
                            (?![a-z])
                        ''', re.X)

    pattern9 = re.compile(r'''                  # "x-9-days" = "for 9 days"
                            (?<![a-z])
                            x                   # mandatory
                            (\s)? (\-\s?)?
                            (?P<time_interval>TIMEINTERVAL)
                            (?P<thereafter>\s
                                (
                                    there \s? after     |
                                    after \s that       |
                                    more
                                )
                            )?
                            (?![a-z])
                        ''', re.X)

    pattern10 = re.compile(r'''                  # "x7" at the end of the sig means "for 7 timeunits" where timeunit needs to be looked up from left context
                            (?<![a-z])
                            (x|for)             # mandatory
                            (\s)? (\-\s?)?
                            (?P<quant_after_x>QUANT)
                            \s? \.? \s? $
                        ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())

        if 'time_interval' in groupname2group:      # pattern0 or pattern6 or pattern9
            time_interval_start = match_object.start('time_interval')
            time_interval = parse.position2struc(time_interval_start)

        elif 'time_unit' in groupname2group:        # pattern1 or pattern2
            time_unit_start = match_object.start('time_unit')
            time_unit = parse.position2struc(time_unit_start)
            quant_start = match_object.start('quant')
            quant = parse.position2struc(quant_start)
            time_interval = TimeInterval(quant, time_unit)

        elif 'time_unit_no_quant' in groupname2group:        # pattern8 - implied quant = 1
            time_unit_start = match_object.start('time_unit_no_quant')
            time_unit = parse.position2struc(time_unit_start)
            quant = Quant(num_type = 'int', value = 1)
            time_interval = TimeInterval(quant, time_unit)

        elif 'nights' in groupname2group:        # pattern3
            if 'plural' in groupname2group:
                plurality = 'plural'
            else:
                plurality = 'singular'
            time_unit = TimeUnit('day', constituents = [match_object.group('nights')], plurality = plurality)
            quant_start = match_object.start('quant')
            quant = parse.position2struc(quant_start)
            time_interval = TimeInterval(quant, time_unit)

        elif 'quant_after_x' in groupname2group:        # pattern10 - implied time unit needs to be looked up from left context
            quant_start = match_object.start('quant_after_x')
            quant = parse.position2struc(quant_start)
            if quant.num_type != 'int':
                return None
            left_pattern = re.compile('(?P<left_struc>PERIODICITY|FREQ)')
            left_obj = left_pattern.search(left_context)
            if not left_obj:
                return None
            left_groupname2group = trim_dictionary(left_obj.groupdict())
            left_struc_start = left_obj.start('left_struc')
            left_struc = parse.position2struc(left_struc_start)
            left_time_unit_value = left_struc.time_unit.value
            if left_struc.label == 'PERIODICITY' and left_time_unit_value == 'hour':      # "1tsp q8h x10"  means "for 10 days'. Duration is probably in hours.
                left_time_unit_value = 'day'
            plurality = 'plural' if quant.value > 1 else 'singular'
            time_unit = TimeUnit(value = left_time_unit_value, constituents = [], plurality = plurality)
            time_unit.rules_used.append('*deduced_from_left_context_in_presence_of_x_Quant')
            time_interval = TimeInterval(quant, time_unit)

        else:           # should not happen
            return None

        if 'terminal_time_interval' in groupname2group:      # pattern5  Make sure it is preceded by a structure that can't be combined with TIMEINTERVAL for anything else
            good_left_context = re.compile('(FREQ|PERIODICITY|TIMING|DOSE|ROUTE|AS_NEEDED|SITE|VEHICLE)\s$')
            if not good_left_context.search(left_context):
                return None

        duration = Duration(time_interval)

        if 'start_next' in groupname2group:
            duration.offset = 'next'
        if 'additional' in groupname2group or 'thereafter' in groupname2group:
            duration.offset = 'next'
        elif 'start_first' in groupname2group:
            start_first_start = match_object.start('start_first')
            start_first = parse.position2struc(start_first_start)
            if start_first.value == 1 and start_first.num_type == 'ord':
                duration.offset = 'first'
        elif 'up_to' in groupname2group:
            duration.up_to = True

        if 'directive_remove' in groupname2group:      # pattern7
            # "Remove after 3 weeks" means "do whatever primary directive says for 3 weeks, then Remove", i.e. turn this to Duration Directive/remove
            # Same with "stop in 3 weeks"
            directive_remove_start = match_object.start('directive_remove')
            directive_remove = parse.position2struc(directive_remove_start)
            if directive_remove.value not in ('remove', 'stop'):
                return None
            else:
                space = Struc(' ', [])
                then = ThenChrono([])
                return [duration, space, then, directive_remove]

        if 'skip' in groupname2group:           # pattern8
            # "skip 7 days" or "skip a week" = "stop for 7 days", etc.
            directive = Directive(value = 'stop')
            directive.constituents = ['skip']
            space = Struc(' ', [])
            return [directive, space, duration]

        return [duration]

    rule = Rule_ExtractStrucs(   name = 'duration',
                                    search_patterns = [pattern0, pattern1, pattern2, pattern3, pattern4, pattern5,
                                                       pattern6, pattern7, pattern8, pattern9, pattern10],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule


def rule_assembly_repeat():
    """
    Repeat instruction: e.g.:
            "may repeat every 5 minutes up to 3 times",
            "take one tablet by mouth at onset of headache may repeat in 2 hours"
            "take 1 tablet at onset of migraine. may repeat once after 2 hours. max 10 mg/day."
            "1 tablet once daily and repeat 1 tablet 1 week later"
            "apply 1 patch once weekly as directed for 3 weeks, leave off for 1 week then repeat cycle"
    """


    pattern0 = re.compile(r'''                              #  Include Periodicity or implied periodicity "5 mins apart" or offset
                            (?<![a-z])
                            (
                                (?P<prev_repeat>REPEAT)                     |
                                ((dosage \s)? may \s be \s)? repeated       |
                                ((you \s)? may \s)? repeat (\s (the \s)? treatment)? (\s (the \s)? cycle)?
                            ) \s
                            (
                                (?P<periodicity>PERIODICITY) \s?                            |
                                (?P<timeinterval_periodicity>TIMEINTERVAL) \s apart \s?     |
                                (
                                    (in \s)? (after \s)?
                                    (?P<offset>TIMEINTERVAL)
                                    (\s later)?
                                )
                            )
                            (?![a-z])
                        ''', re.X)

    pattern1 = re.compile(r'''                  # Include reps or max_reps ("repeat up to 3 times" or "repeat 3 times" "or repeate for 3 doses"
                            (?<![a-z])
                            (
                                (?P<prev_repeat>REPEAT)                     |
                                (?P<periodicity>PERIODICITY)                |   # even without word "repeat", "take 1 tab every 5 mins for 3 doses
                                ((dosage \s)? may \s be \s)? repeated       |
                                ((you \s)? may \s)? repeat (\s (the \s)? treatment)? (\s (the \s)? cycle)?
                            ) \s
                            (?P<reps_or_max_reps>
                                (                               # optional
                                    (for \s)?
                                    (?P<upto>
                                        up (\s|\-)? to      |
                                        max (\s of)?
                                    ) \s?                                   |
                                    for (\s a \s total (\s of)?)? \s?
                                )?
                                (?P<reps>QUANT)                     # mandatory
                                \/?\s?
                                (more \s)?
                                (times? | doses?)           # mandatory
                            )
                            (?![a-z])
                        ''', re.X)

    pattern2 = re.compile(r'''                  # "up to 3 doses", even without "repeat"
                            (?<![a-z])
                            (
                                (
                                    (?P<prev_repeat>REPEAT)                     |
                                    ((dosage \s)? may \s be \s)? repeated       |
                                    ((you \s)? may \s)? repeat (\s (the \s)? treatment)? (\s (the \s)? cycle)?
                                ) \s
                            )?
                            (                                   # mandatory
                                (for \s)?
                                (?P<upto>
                                    up (\s|\-)? to  |
                                    max (\s of)?
                                )                           |
                                for (\s a)? \s total \s of
                            )
                            \s?
                            (?P<reps>QUANT)
                            \/?\s?
                            (more \s)?
                            doses                           # the word "doses" is mandatory if there is no clear "repeat"
                            (?![a-z])
                        ''', re.X)

    pattern3 = re.compile(r'''                  # Include condition before REPEAT, eg "if no relief, repeat ..."
                            (?<![a-z])
                            (?P<condition>
                                (?P<no_relief>
                                    (if \s (head \s? ache \s)? )?
                                    (no \s relief | not \s resolved | not \s improve(s|d) (\s or \s recurrs)?| condition \s persists?)
                                )
                            ),? \s
                            (
                                (?P<prev_repeat>REPEAT)                     |
                                ((dosage \s)? may \s be \s)? repeated       |
                                ((you \s)? may \s)? repeat (\s (the \s)? treatment)? (\s (the \s)? cycle)?
                            )
                            (?![a-z])
                        ''', re.X)

    pattern4 = re.compile(r'''                  # Include condition after REPEAT eg "repeat up to 5 times if condition persists"
                            (?<![a-z])
                            (
                                (?P<prev_repeat>REPEAT)                     |
                                ((dosage \s)? may \s be \s)? repeated       |
                                ((you \s)? may \s)? repeat (\s (the \s)? treatment)? (\s (the \s)? cycle)?
                            ) \s
                            (?P<condition>
                                (?P<no_relief>
                                    (if \s (head \s? ache \s)? )?
                                    (no \s relief | not \s resolved | not \s improve(s|d) (\s or \s recurrs)?| condition \s persists?)
                                )
                            )
                            (?!\sMISCELLANEOUS)     # the condition may belong to Call_911, in which case omit it here (it will be repeated in MISC typ = call_911
                            (?!\scall\sQUANT)       # the condition may belong to Call_911, in which case omit it here (it will be repeated in MISC typ = call_911
                            (?![a-z])
                        ''', re.X)

    pattern5 = re.compile(r'''                  # Then Repeat = Repeat, And repeat = Repeat
                            (THEN_CHRONO |  AND_CONJ)
                            \s
                            (
                                (?P<prev_repeat>REPEAT)                     |
                                ((dosage \s)? may \s be \s)? repeated       |
                                (start \s over)                             |
                                ((you \s)? may \s)? repeat (\s (the \s)? treatment)? (\s (the \s)? cycle)?
                            )
                            (?![a-z])
                        ''', re.X)



    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())
        periodicity = None
        offset = None
        reps = None
        max_reps = None
        condition = None
        if 'periodicity' in groupname2group:
            periodicity_start = match_object.start('periodicity')
            periodicity = parse.position2struc(periodicity_start)
        elif 'timeinterval_periodicity' in groupname2group:
            timeinterval_periodicity_start = match_object.start('timeinterval_periodicity')
            timeinterval_periodicity = parse.position2struc(timeinterval_periodicity_start)
            quant = timeinterval_periodicity.quant
            time_unit = timeinterval_periodicity.time_unit
            periodicity = Periodicity(quant, time_unit)
        elif 'offset' in groupname2group:
            offset_start = match_object.start('offset')
            offset = parse.position2struc(offset_start)


        if 'reps' in groupname2group:
            reps_start = match_object.start('reps')
            reps = parse.position2struc(reps_start)
            if reps.num_type not in ('int', 'var'):
                # E.g. "take 1 tab daily for first dose" should be excluded.
                return None
            if 'upto' in groupname2group:
                max_reps = reps
                reps = None

        if 'condition' in groupname2group:
            if 'no_relief' in groupname2group:
                condition = 'no_relief'

        if 'prev_repeat' in groupname2group:
            prev_repeat_start = match_object.start('prev_repeat')
            prev_repeat = parse.position2struc(prev_repeat_start)
            periodicity = periodicity if periodicity else prev_repeat.periodicity
            offset = offset if offset else prev_repeat.offset
            reps = reps if reps else prev_repeat.reps
            max_reps = max_reps if max_reps else prev_repeat.max_reps
            condition = condition if condition else prev_repeat.condition
            repeat = Repeat(periodicity = periodicity, offset = offset, reps = reps, max_reps = max_reps, condition = condition)
            return [repeat]

        repeat = Repeat(periodicity, offset, reps, max_reps, condition)
        repeat.constituents = [match_object.group()]

        return [repeat]

    rule = Rule_ExtractStrucs(   name = 'repeat',
                                    search_patterns = [pattern0, pattern1, pattern2, pattern3, pattern4, pattern5],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule

def rule_assembly_taper():
    """
    Taper instruction: e.g.:
            "take one capsule three times daily may increase by 1 cap every 3 days to effective dose or max three three times daily",
            "start with 1 capsule at bedtime then increase by 1 capsule a day to effect"
            "inject 10 units subcutaneously every morning and increase by 2 units every other day until fbs 80-120"
            "take one tablet every 12 hours -increase to 2 tablets every 12 hours after 1 week "
            "take 1 tablet by mouth each night at bedtime may increase to 2 tablet"
            "take 6 tablets by mouth on day 1 then decrease by 1 tablet every 3 days"
    """


    pattern0 = re.compile(r'''                              #   Basic pattern: capture typ, direction, dose: e.g. "increase by 1 tab"
                            (?<![a-z])
                                ( (may|can) \s)?
                                (?P<direction>          # Mandatory: increase vs dicrease.
                                    (?P<increase>
                                        increase
                                    )                   |
                                    (?P<decrease>
                                        decrease
                                    )
                                )
                                (\s (the \s)? dose)?
                                (                                   # optional (primarily to parse the dictionary items which is "MAY INCREASE THE DOSE EVERY DAY BY <<NUM_0>> CAPSULE(S)"
                                    \s (?P<increment_periodicity>PERIODICITY)
                                )?
                                \s
                                (?P<taper_type>         # Optional: increase/decrease BY this amount every 2 days vs. increase/decrease TO this amount after 1 week
                                                        # If taper_type is absent, assume incremental. E.g. "increase 1 capsule" in the data seems to imply "increase by 1 capsule"
                                    (
                                        (?P<incremental>
                                            by
                                        )                   |
                                        (?P<target>
                                            to
                                            (\s (?P<optional_directive>DIRECTIVE))?     # frequent expression: "may increase to TAKE 2 tabs ..."
                                        )
                                    )
                                    \s
                                )?
                                (?P<dose_or_quant>                       # Either dose or quantity (from which we will deduce dose) is mandatory
                                    (?P<dose>DOSE)          |
                                    (?P<quant_dose>QUANT)
                                )
                            (?![a-z])
                        ''', re.X)

    pattern1 = re.compile(r'''                              #   Capture dincrement_periodicity
                            (?<![a-z])
                                (?P<prev_taper>                             # Mandatory
                                    TAPER
                                )
                                (\s (?P<route>ROUTE))?      # Optional this is usually junk -- we should not be telling the patient for the first time in the tapering instructions what the route should be.
                                \s
                                (?P<periodicity>PERIODICITY)    # Mandatory
                            (?![a-z])
                        ''', re.X)

    pattern2 = re.compile(r'''                              #   Capture stop_condition (e.g. "to effect"
                            (?<![a-z])
                                (?P<prev_taper>          # Mandatory
                                    TAPER
                                )
                                \s
                                (?P<stop_condition>      # Mandatory
                                    STOP_CONDITION
                                )
                            (?![a-z])
                        ''', re.X)


    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())

        direction = None
        taper_type = None
        dose = None
        dose_freq = None
        dose_periodicity = None
        increment_periodicity = None
        offset = None
        stop_condition = None

        if 'prev_taper' in groupname2group:
            prev_taper_start = match_object.start('prev_taper')
            prev_taper = parse.position2struc(prev_taper_start)
        else:
            prev_taper = None

        if 'direction' in groupname2group:
            if 'increase' in groupname2group:
                direction = 'increase'
            else:
                direction = 'decrease'

        if 'taper_type' in groupname2group:
            if 'incremental' in groupname2group:
                taper_type = 'incremental'
            else:
                taper_type = 'target'
        elif not prev_taper:
            # If taper_type is absent, assume incremental. E.g. "increase 1 capsule" in the data seems to imply "increase by 1 capsule"
            taper_type = 'incremental'

        if 'dose_or_quant' in groupname2group:
            if 'dose' in groupname2group:
                dose_start = match_object.start('dose')
                dose = parse.position2struc(dose_start)
            else:
                # it's just a number, e.g. "increase by 1 every day'
                quant_dose_start = match_object.start('quant_dose')
                quant_dose = parse.position2struc(quant_dose_start)
                left_pattern = re.compile('(?P<left_dose>DOSE)')
                left_dose = left_pattern.search(left_context)
                if left_dose:
                    left_dose_start = left_dose.start('left_dose')
                    left_dose = parse.position2struc(left_dose_start)
                    form_value = left_dose.form.value
                    if quant_dose.value > 1:
                        plurality = 'plural'
                    else:
                        plurality = 'singular'
                    form = Form(form_value, plurality, [])
                    form.rules_used.append('*deduced_in_taper_identification*')
                else:
                    form = None
                dose = Dose(quant = quant_dose, form = form)

        if 'increment_periodicity' in groupname2group:
            increment_periodicity_start = match_object.start('increment_periodicity')
            increment_periodicity = parse.position2struc(increment_periodicity_start)

        if 'periodicity' in groupname2group:
            periodicity_start = match_object.start('periodicity')
            periodicity = parse.position2struc(periodicity_start)
            if prev_taper and prev_taper.taper_type == 'target':
                dose_periodicity = periodicity
            else:
                increment_periodicity = periodicity

        if 'stop_condition' in groupname2group:
            stop_condition_start = match_object.start('stop_condition')
            stop_condition = parse.position2struc(stop_condition_start)

        if prev_taper:
            if dose_freq:
                prev_taper.dose_freq = dose_freq
            if dose_periodicity:
                prev_taper.dose_periodicity = dose_periodicity
            if increment_periodicity:
                prev_taper.increment_periodicity = increment_periodicity
            if offset:
                prev_taper.offset = offset
            if stop_condition:
                prev_taper.stop_condition = stop_condition
            return [prev_taper]

        if 'THEN_CHRONO' in left_context or 'then' in left_context:
            then_flag = True
        else:
            then_flag = False

        taper = Taper(taper_type = taper_type, direction = direction, dose = dose)
        taper.dose_freq = dose_freq
        taper.dose_periodicity = dose_periodicity
        taper.increment_periodicity = increment_periodicity
        taper.offset = offset
        taper.stop_condition = stop_condition
        taper.then_flag = then_flag
        taper.constituents = [match_object.group()]

        return [taper]

    rule = Rule_ExtractStrucs(   name = 'taper',
                                    search_patterns = [pattern0, pattern1, pattern2],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule




def rule_assembly_timing():

    pattern0 = re.compile(r'''
                                (?<![a-z])

                                (?P<start_flag>                         # Optional "start before the procedure". Can be combined with offset. e.g. "to start 3 days prior to surgery",
                                    (to \s start \s)    |               # or work on it's own: "start after surgery"
                                    (to \s begin \s)    |
                                    (start(ing)? \s)    |
                                    (begin(ning)? \s)
                                )?
                                (                                       # optional offset e.g. "2 hours before lunch" or "within 30 minutes of onset of severe migraine"
                                    (within \s)?
                                    (?P<offset>TIMEINTERVAL) \s
                                )?
                                ((?P<relation>                              # mandatory if the left context doesn't have "TIMING and" immediately preceding
                                            in                          |
                                            at (\s (each|every))?       |
                                            before (\s (each|every|at))?        |    # there is a bizarre use "take 2 hours before at bedtime" ("before at")
                                            prior (\s to)? (\s (each|every))?   |
                                            after (\s (each|every))?    |
                                            during (\s (each|every))?   |
                                            every                       |           # E.g. "take with each meal"
                                            each                        |
                                            on                          |
                                            of                          |   # e.g. "within 30 mins of onset of headache
                                            upon                        |
                                            with (\s (each|every))?
                                ) \s)?
                                ((the|an|a) \s)?
                                (?P<landmark>                               # There are 3 types of landmarks: meal events, day parts, and events.
                                        (?P<breakfast>
                                            breakfast               |
                                            morning \s meals?
                                        )                           |
                                        (?P<lunch>
                                            lunch
                                        )                           |
                                        (?P<dinner>
                                            dinner                  |
                                            supper                  |
                                            evening \s meals?       |
                                            last \s meal \s of \s the \s day
                                        )                           |
                                        (?P<empty_stomach>
                                            empty \s
                                            (stomach|stoamch)
                                        )                           |
                                        (?P<food_or_milk>
                                            food \s or \s milk
                                        )                           |
                                        (?P<food>
                                            food    |
                                            eating  |
                                            you \s eat
                                        )                           |
                                        (?P<meals>meals?)           |       # End of meal events
                                        (?P<evening>evenings?)      |           # Start of day parts
                                        (?P<pm>
                                            p\.?m\.?
                                        )                           |
                                        (?P<morning>mornings?)      |
                                        (?P<am>a\.?m\.?)            |
                                        (?P<afternoon>after \s? noons?) |
                                        (?P<bedtime>
                                            bed(\s|\-)?time                 |
                                            going \s to \s (bed|sleep)      |
                                            bed                             |
                                            sleep
                                        )                           |
                                        (?P<awakening>a?wakening)   |
                                        (?P<night>
                                            nights?     |
                                            nite
                                        )                           |           # End of day parts
                                        (?P<sex>                                # Start of events
                                            (having \s)?
                                            sex (ual)?
                                            (\s activity)?
                                            (\s intercourse)?
                                            (\s relations?)?
                                        )                           |
                                        (?P<onset_of_headache>
                                            (
                                                onset \s of \s (severe \s)?
                                                    (
                                                        head \s ? ache      |
                                                        migraine?s?         |
                                                        (?P<indication_headache>INDICATION)
                                                    )
                                                |
                                                (severe \s)? (migraine?s? | head \s? aches?) \s onset
                                            )
                                        )                           |
                                        (?P<bowel_movement>             # after each bowel movement
                                            (?<!for\s)
                                            (?<!loose\s)
                                            bowels? \s movements?
                                        )                           |
                                        (?P<loose_stool>                # "after each loose stool" or "after loose bowel movement"
                                            (?<!for\s)
                                            (
                                                loose   |
                                                lose
                                            )
                                            \s
                                            (
                                                stool   |
                                                bowels? \s movements?
                                            )
                                        )                           |
                                        (?P<diaper_change>
                                            diaper \s
                                            change
                                        )                           |
                                        (?P<dental_appt>
                                            dental \s
                                            (appointment|appt|visit|exam(ination)?|work)
                                        )                           |
                                        (?P<appt>
                                            (your \s)?
                                            (
                                                appointment     |
                                                appt            |
                                                app             |
                                                appointmant
                                            )
                                        )                           |
                                        (?P<procedure>
                                            (dental \s)?
                                            (
                                                procedures?     |
                                                surgery         |
                                                treatments?     |
                                                biopsy          |
                                                colonoscopy     |
                                                dialysis
                                            )
                                        )
                                )
                                (?![a-z])
                            ''', re.X)

    pattern1 = re.compile(r'''      # special cases
                                (?P<special_cases>
                                    (?P<awakening>when \s you \s wake \s up)    |
                                    (?P<day_prior_surgery>                      # To catch "10 pm day prior to surgery"
                                        (?P<time_of_day>TIMING)
                                        \s
                                        ((on \s)? the \s)?
                                        (?P<timeunit_day>TIMEUNIT)
                                        \s
                                        (?P<procedure_timing>TIMING)
                                    )
                                )

                            ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        """
        Processes timing indicators such as "1 hour before every meal", as well as mulitple indicators ("after breakfast and dinner").

        Multiple indicators has to be signaled by "TIMING and" in the leftmost context. This is the only case where relation is not mandatory.

        If the landmark is a daypart (evening, morning, afternoon, bedtime, night) the relation doesn't matter, it is purely
        a linguistic convenience. We omit the relation then (relation = None).

        If the landmark is a meal type (breakfast, lunch, dinner, food, meal) then there are 3 types of semantic relations: before, during, after.

        "every_flag" is a T/F flag representing that the action should take place at the specified landmark at each of it's occurences.
        E.g. "with each meal and at bedtime", "every night at bedtime".
        If the landmark is a once-a-day event (breakfast, bedtime, etc) then the presence of "every" or "each" implies that
        the action should happen every day, i.e. we need to add something like "Take this medicine every day" if
        it is not expressed elsewhere in the Sig. But if the landmark can happen several times daily, e.g. "after every meal",
        then the daily freq is possibly different from 1.
        """

        groupname2group = trim_dictionary(match_object.groupdict())

        if 'special_cases' in groupname2group:      # pattern1
            if 'awakening' in groupname2group:
                landmark = 'awakening'
                timing = Timing(landmark = landmark, relation = 'at', offset = None, start_flag = False, every_flag = False)
                return [timing]
            elif 'day_prior_surgery' in groupname2group:
                time_of_day_start = match_object.start('time_of_day')
                time_of_day = parse.position2struc(time_of_day_start)
                timeunit_day_start = match_object.start('timeunit_day')
                timeunit_day = parse.position2struc(timeunit_day_start)
                procedure_timing_start = match_object.start('procedure_timing')
                procedure_timing = parse.position2struc(procedure_timing_start)
                if time_of_day.typ != 'time_of_day' or timeunit_day.value != 'day' or procedure_timing.offset or procedure_timing.typ != 'event':
                    return
                quant = Quant(num_type = 'int', value = 1)
                procedure_timing.offset = TimeInterval(quant = quant, time_unit = timeunit_day)
                space = Struc(label = ' ')
                return [time_of_day, space, procedure_timing]
            return None

        left_context_timing_found = re.search("(?P<prior_timing>TIMING) (and|AND_CONJ)\s?$", left_context)
        if not left_context_timing_found and not 'relation' in groupname2group:
            return None
        elif left_context_timing_found and not 'relation' in groupname2group:
            # In a situation like "take 2 tabs with breakfast and dinner", when analyzing "dinner" use the relation from "breakfast".
            prior_timing_start = left_context_timing_found.start('prior_timing')
            prior_timing = parse.position2struc(prior_timing_start)
            relation = prior_timing.relation if prior_timing.relation else ''
        else:
            relation = groupname2group['relation']

        if 'every' in relation or 'each' in relation:       # e.g. each meal
            every_flag = True
        else:
            every_flag = False

        relation = relation.replace('each', '').replace('every', '')
        relation = trim(relation)

        if relation in ('before', 'prior', 'prior to'):
            relation = 'before'
        elif relation in ('at', 'during', 'with', 'upon', 'of'):
            relation = 'at'
        elif relation == 'after':
            relation = 'after'
        else:
            relation = ''


        landmark = groupname2group['landmark']
        landmark = landmark.replace('.', '')

        if 'morning' in groupname2group and 'and evening meal' in right_context:
            # there is a frequent expression morning and evening meals
            landmark = 'breakfast'
        elif 'am' in groupname2group:
            landmark = 'morning'
        elif 'afternoon' in groupname2group:
            landmark = 'afternoon'
        elif 'appt' in groupname2group:
            landmark = 'appointment'
        elif 'awakening' in groupname2group:
            landmark = 'awakening'
        elif 'bedtime' in groupname2group:
            landmark = 'bedtime'
        elif 'bowel_movement' in groupname2group:
            landmark = 'bowel_movement'
        elif 'breakfast' in groupname2group:
            landmark = 'breakfast'
        elif 'dental_appt' in groupname2group:
            landmark = 'dental appointment'
        elif 'diaper_change' in groupname2group:
            landmark = 'diaper change'
        elif 'dinner' in groupname2group:
            landmark = 'dinner'
        elif 'empty_stomach' in groupname2group:
            landmark = 'empty stomach'
        elif 'evening' in groupname2group:
            landmark = 'evening'
        elif 'food_or_milk' in groupname2group:
            landmark = 'food or milk'
        elif 'food' in groupname2group:
            landmark = 'food'
        elif 'loose_stool' in groupname2group:
            landmark = 'loose_stool'
        elif 'lunch' in groupname2group:
            landmark = 'lunch'
        elif 'meals' in groupname2group:
            landmark = 'meal'
        elif 'morning' in groupname2group:
            landmark = 'morning'
        elif 'night' in groupname2group:
            landmark = 'night'
        elif 'procedure' in groupname2group:
            if not relation:            # E.g. we can have expressions such as "instill in surgery eye". To be Timing for procedure or surgery, need to have a relation.
                return None
            landmark = 'procedure'
        elif 'onset_of_headache' in groupname2group:
            if 'indication_headache' in groupname2group:
                indication_start = match_object.start('indication_headache')
                indication = parse.position2struc(indication_start)
                if indication.condition not in ('headache', 'migraine'):
                    return None
            landmark = 'onset_of_headache'
        elif 'pm' in groupname2group:
            landmark = 'evening'
        elif 'sex' in groupname2group:
            landmark = 'sex'


        if 'offset' in groupname2group:
            offset_start = match_object.start('offset')
            offset = parse.position2struc(offset_start)
        else:
            offset = None

        if 'start_flag' in groupname2group:
            start_flag = True
        else:
            start_flag = False


        timing = Timing(landmark = landmark, relation = relation, offset = offset, start_flag = start_flag, every_flag = every_flag)

        # "every monring" means "every day in the morning", so include Periodicity of 1 day
        if every_flag and landmark in ('morning', 'evening', 'night', 'afternoon', 'bedtime', 'lunch'):
            time_unit = TimeUnit(value = 'day', constituents = [], plurality = 'singular')
            quant = Quant(num_type = 'int', value = 1)
            periodicity = Periodicity(quant = quant, time_unit = time_unit)
            space = Struc(label = ' ')
            return [periodicity, space, timing]

        return [timing]

    rule = Rule_ExtractStrucs(   name = 'timing',
                                    search_patterns = [pattern0, pattern1],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule


def rule_assembly_time_of_day():

    pattern0 = re.compile(r'''     # searches for time of day, e.g. "at 12:00 noon", "at five p.m.", or "at noon".  "at" is mandatory, o/w we will wrongly process Latin sigs such as "t1amqd"
                                (?<![a-z])
                                (                               # Mandatory: "at" or "around" or "every" parenthesis is mandatory. E.g. "at 8 am" or "(8 am, 10pm)"
                                    (at \s)? (?P<every> every \s)   |        # e.g "at every 6pm" which implies "every day at 6pm"
                                    at \s                           |
                                    @ \s?                           |
                                    \(                              |
                                    around \s
                                )
                                (?P<time_of_day>
                                    (
                                        (?P<noon>noon)                          |
                                        (
                                            (
                                                (?P<hour_alone>QUANT)           |
                                                (?P<hour>QUANT):(?P<min>QUANT)
                                            )
                                            \s?
                                            (
                                                (?P<am>a\.?m\.?)                |
                                                (?P<pm>p\.?m\.?|noon)           |   # for expression such as at "12 noon"
                                                (o\'?clock)
                                            )
                                        )
                                    )
                                )
                                (?![a-z])
                            ''', re.X)

    pattern1 = re.compile(r'''          # Without "at" upfront, but then should have comma at the end and space upfront.
                                        # E.g. "take one tablet by mouth 9 am", "take one tablet by mouth 1pm, 2pm, and 10pm day prior to surgery"
                                        # Still, because of possible confusion with latin "1am" or "1pm" meaning "take 1 in the morning" or "take 1 in the evening", we make sure that
                                        # Quant is greater than 3 (doses are usually < 4 while times are usually in the morning >6am or in the evening after 6pm)
                                (?<=\s)
                                (?P<time_of_day_sans_at>
                                    (
                                        (?P<hour_alone>QUANT)           |
                                        (?P<hour>QUANT):(?P<min>QUANT)
                                    )
                                    \s?
                                    (
                                        (?P<am>a\.?m\.?)                |
                                        (?P<pm>p\.?m\.?|noon)               # for expression such as at "12 noon"
                                    )
                                )
                                (,|$)                                       # has to have comma or be at the end
                            ''', re.X)

    pattern2 = re.compile(r'''      # searches for time of day in a sequence e.g. "at 9 am, 12 noon, and 5 pm" where there is no "at" before QUANT
                                    # In other words, if we have "TIMING, and 5 pm. If "and" is missing, need to add it.
                                (?<=TIMING)
                                (?P<time_of_day_in_sequence> (\s|,) \s? (AND_CONJ)? \s?)        # regardless of whether AND_CONJ is there, we will need to (re)insert one.
                                (?P<time_of_day>
                                    (
                                        (?P<noon>noon)                          |
                                        (
                                            (
                                                (?P<hour_alone>QUANT)           |
                                                (?P<hour>QUANT):(?P<min>QUANT)
                                            )
                                            \s?
                                            (
                                                (?P<am>a\.?m\.?)                  |
                                                (?P<pm>p\.?m\.?|noon)                # for expression such as at "12 noon"
                                            )
                                        )

                                    )

                                )
                                (?![a-z])
                            ''', re.X)

    pattern3 = re.compile(r'''      # For processing Dictionary item only (with period at the end)
                        (?<=\s)
                        at \s (?P<dict_time>QUANT) (\s o\'?clock)\.
                        ''', re.X)


    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        """
        """

        groupname2group = trim_dictionary(match_object.groupdict())
        keys = groupname2group.keys()

        if 'noon' in groupname2group:
            landmark = '12:00'
        elif 'dict_time' in groupname2group:    # pattern3 only
            quant_start = match_object.start('dict_time')
            quant = parse.position2struc(quant_start)
            if quant.num_type == 'var' and quant.var_type == 'time':
                landmark = '<<' + quant.value_repr + '>>'
            else:
                return
        else:
            if 'am' in groupname2group:
                am = True
            else:
                am = False
            if 'hour_alone' in groupname2group:
                hour_start = match_object.start('hour_alone')
                hour = parse.position2struc(hour_start)
                minutes = '00'
            elif 'hour' in groupname2group:
                hour_start = match_object.start('hour')
                hour = parse.position2struc(hour_start)
                min_start = match_object.start('min')
                minutes = parse.position2struc(min_start)
                if minutes.num_type != 'int' or minutes.value >=60:
                    return None
                minutes = ('0' + str(minutes.value))[-2:]
            if hour.num_type != 'int':
                return None
            hour = hour.value
            if am and hour > 12:
                return None
            elif hour > 24:
                return None
            if not am and hour < 12:            # use military time, so if pm and hour < 12, add 12.
                hour += 12
            hour = ('0' + str(hour))[-2:]
            landmark = hour + ':' + minutes
        relation = ''
        offset = None
        every_flag = None
        typ = 'time_of_day'
        timing = Timing(landmark, relation, typ, offset, every_flag)

        if 'time_of_day_in_sequence' in groupname2group:    # pattern2 only
            space1 = Struc(label = ' ')
            space2 = Struc(label = ' ')
            and_conj = AndConj(constituents = [''])
            return [space1, and_conj, space2, timing]

        if 'every' in groupname2group:    # pattern0 and pattern1
            quant = Quant(num_type = 'int', value = 1)
            every_day = Periodicity(quant, TimeUnit('day', []))
            space = Struc(label = ' ')
            return [every_day, space, timing]

        if 'time_of_day_sans_at' in groupname2group and 'hour_alone' in groupname2group and not right_context:           # pattern1
            # To avoid possible confusion with latin "1am" or "1pm" meaning "take 1 in the morning" or "take 1 in the evening", we make sure that
            # Quant is greater than 3 (doses are usually < 4 while times are usually in the morning >6am or in the evening after 6pm)
            # if there is no comma at the end.
            hour_int = int(hour)
            if hour_int < 4:
                timing.landmark = 'morning'
                timing.typ = 'day_part'
            elif hour_int > 12 and hour_int < 16:
                timing.landmark = 'evening'
                timing.typ = 'day_part'


        return [timing]

    rule = Rule_ExtractStrucs(   name = 'timing_time_of_day',
                                    search_patterns = [pattern0, pattern1, pattern2, pattern3],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule



def rule_assembly_calendar_event():
    """ Finds Calendar_Event strucs. They are of the type
            Now, today, the first day, 1st day, Day 1, days 2-5, monday

    """

    pattern0 = re.compile(r'''      # searches for "days 1 through 5" or "on week 1" or "day-1 to day-7"
                                    # But if there is a dash between TimeUnit and Quant, make sure it is not part of "max per day-8"
                                (?<![a-z])
                                (?P<prefix> (on|for) \s)?       # optional
                                (?P<timeunit>TIMEUNIT) (\s|\-)? (?P<quant>QUANT)
                                (?!\s\-|\-)
                            ''', re.X)

    pattern1 = re.compile(r'''     # searches for "the first day"
                                (?<![a-z])
                                ((on \s the|on|the) \s)?    # could start with "on the" or "on" or "the", but may not
                                                            # E.g. take 2 capsules every 12 hours first day, then 1 ..."
                                (?P<quant_ord>QUANT) \s? (?P<timeunit>TIMEUNIT)
                            ''', re.X)


    pattern2 = re.compile(r'''     # searches for now, at once (which we change to "at 1 time"), today, tomorrow
                                (?<![a-z])
                                (?P<relative>
                                    (?P<now>now)            |   # was (?P<now>now|(?<!all \s)at \s once|(?<!all \s) at \s (?P<quant>QUANT) \s time)
                                                                # but decided that "at once" means as_single_dose (which is a DIR_ATTRIBUTE)
                                    (?P<today>today)        |
                                    (?P<tomorrow>tomorrow)
                                )
                                (?![a-z])
                            ''', re.X)


    pattern3 = re.compile(r'''      # Exceptional case: searches for "the first day" but where "day" is omitted/implied
                                    # e.g. "take 2 tablets by mouth on the first, then take 1 tablet daily until gone")
                                (?<![a-z])
                                (?P<omitted_day>
                                    on \s the \s (?P<quant_ord>QUANT) ,?
                                )
                                (?!\sTIMEUNIT)
                                (?!\sday)
                                (?![a-z])
                            ''', re.X)


    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        """
        """

        groupname2group = trim_dictionary(match_object.groupdict())

        if 'prefix' in match_object.groupdict() and 'prefix' not in groupname2group:    # pattern0. Make sure you are not within the MaxDose framework (e.g. "max per day-8")
            bad_left_context = re.compile(r'per\s$')
            if bad_left_context.search(left_context):
                return None

        if 'timeunit' in groupname2group:               # it's in the form ("on day 1" or "on days 1-5") or "first day"
            time_unit_start = match_object.start('timeunit')
            time_unit = parse.position2struc(time_unit_start)
            # it then has to have either 'quant' or 'quant_ord'
            if 'quant' in groupname2group:
                quant_start = match_object.start('quant')
                quant = parse.position2struc(quant_start)
            elif 'quant_ord' in groupname2group:        # of the form "first day" or "2nd day"
                quant_start = match_object.start('quant_ord')
                quant = parse.position2struc(quant_start)
                if quant.num_type != 'ord':     # this should be a TIMEINTERVAL, not a calendar_event.
                    return None
            calendar_event = Calendar_Event(typ = 'numeric')
            calendar_event.quant = quant
            calendar_event.time_unit = time_unit
            calendar_event.constituents = [time_unit, quant]
        elif 'relative' in groupname2group:
            calendar_event = Calendar_Event(typ = 'relative')
            if 'now' in groupname2group:
                if 'quant' in groupname2group:                  # we have "at one time", so Quant better be 1 and integer.
                    quant_start = match_object.start('quant')
                    quant = parse.position2struc(quant_start)
                    if quant.value != 1 or quant.num_type != 'int':
                        return None
                calendar_event.value = 'now'
            elif 'today' in groupname2group:
                calendar_event.value = 'today'
            elif 'tomorrow' in groupname2group:
                calendar_event.value = 'tomorrow'
            calendar_event.constituents = [match_object.group()]

        elif 'omitted_day' in groupname2group:          # pattern3:  "the first day" but where "day" is omitted/implied
                                                        # e.g. "take 2 tablets by mouth on the first, then take 1 tablet daily until gone")
            quant_start = match_object.start('quant_ord')
            quant = parse.position2struc(quant_start)
            if quant.num_type != 'ord':
                return None
            correct_right_context_pattern = re.compile(r'^,?\s?(THEN_CHRONO|then)')
            if not correct_right_context_pattern.search(right_context):
                return None
            time_unit = TimeUnit(value = 'day', constituents = [], plurality = 'singular')
            calendar_event = Calendar_Event(typ = 'numeric')
            calendar_event.quant = quant
            calendar_event.time_unit = time_unit
            calendar_event.constituents = [time_unit, quant]
        else:
            return None


        return [calendar_event]

    rule = Rule_ExtractStrucs(   name = 'calendar_event',
                                    search_patterns = [pattern0, pattern1, pattern2, pattern3],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule



def rule_assembly_specific_day():
    """ Finds Specific_Day strucs. They are of the type monday, weekend, 1/2/2013


    """

    pattern0 = re.compile(r'''     # searches for day names
                                (?<![a-z])
                                (on \s)?
                                (?P<day_of_week>
                                    (?P<monday>
                                        mon(days?)?
                                    )                           |
                                    (?P<tuesday>
                                        tue(sdays?)?
                                    )                           |
                                    (?P<wednesday>
                                        wed(nesdays?)?
                                    )                           |
                                    (?P<thursday>
                                        thursdays?  |
                                        thurs       |
                                        thur?
                                    )                           |
                                    (?P<friday>
                                        fri(days?)?
                                    )                           |
                                    (?P<saturday>
                                        sat(urdays?)?
                                    )                           |
                                    (?P<sunday>
                                        sun(days?)?
                                    )                           |
                                    (?P<weekdays>
                                        weekdays?
                                    )                           |
                                    (?P<weekend>
                                        weekends?
                                    )
                                )
                                (                           # optional day_part : e.g. "every sunday morning"
                                    \s
                                    (?P<day_part>
                                        morning | evening | afternoon | night | bedtime
                                    )
                                )?
                                (?![a-z])
                            ''', re.X)


    pattern1 = re.compile(r'''     # searches for day name ranges: Monday - Friday, mon thru fri
                                (?P<start_range>SPECIFIC_DAY)
                                \s?
                                (                           # must have to_range (e.g. mon to fri) or or_range (sat or sunday)
                                    (?P<to_range>
                                        (to|thru|through|\-)
                                    )                           |
                                    (?P<or_range>
                                        or (\s on)?
                                    )                           |
                                    (?P<and_range>
                                        (and|AND_CONJ)
                                    )
                                )
                                \s?
                                (?P<end_range>SPECIFIC_DAY)
                                (?![a-z])
                            ''', re.X)


    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        """
        """

        groupname2group = trim_dictionary(match_object.groupdict())
        if 'day_of_week' in groupname2group:
            keys = set(groupname2group.keys())
            keys = list(keys - set(['day_of_week', 'day_part']))
            day_of_week = keys[0]
            specific_day = Specific_Day(typ = 'day_of_week', day_of_week = day_of_week)
        elif 'start_range' in groupname2group:
            start_of_start_range = match_object.start('start_range')
            start_range = parse.position2struc(start_of_start_range)
            start_of_end_range = match_object.start('end_range')
            end_range = parse.position2struc(start_of_end_range)
            if not start_range.day_of_week or not end_range.day_of_week:
                return None
            if start_range.day_of_week == 'monday' and end_range.day_of_week == 'friday' and 'to_range' in groupname2group:
                specific_day = Specific_Day(typ = 'day_of_week', day_of_week = 'weekdays')
            elif start_range.day_of_week == 'monday' and end_range.day_of_week == 'saturday' and 'to_range' in groupname2group:
                specific_day_0 = Specific_Day(typ = 'day_of_week', day_of_week = 'weekdays')
                specific_day_1 = Specific_Day(typ = 'day_of_week', day_of_week = 'saturday')
                # put commas to prevent infinite loop
                return [specific_day_0, Struc(','), AndConj([]), Struc(','), specific_day_1]
            elif start_range.day_of_week == 'saturday' and end_range.day_of_week == 'sunday':
                specific_day = Specific_Day(typ = 'day_of_week', day_of_week = 'weekend')
            elif 'or_range' in groupname2group:
                # for things like "sunday or monday" just choose one.
                specific_day = start_range
            else:
                return None
        else:
            return None

        if 'day_part' in groupname2group:
            #  e.g. "every sunday morning"
            day_part = groupname2group['day_part']
            timing = Timing(landmark = day_part, relation = '', offset = None, start_flag = False, every_flag = False)
            return [specific_day, Struc(' ', []), timing]

        return [specific_day]

    rule = Rule_ExtractStrucs(   name = 'specific_day',
                                    search_patterns = [pattern0, pattern1],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule



def rule_assembly_nightly():
    """
    REMOVED 11/29/2012
    Now this function is performed within the Periodicity rule.
    """

    pattern = re.compile(r'''
                            (?<![a-z])
                            (?P<nightly>
                                (every|each) \s night   |
                                nightly
                            )
                            (
                                \s(?P<timing>TIMING)    |
                                \s?\.?$
                            )
                        ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        """
        The phrase is either "nightly at bedtime" or "... nightly."
        In the first case, it just means "every day at bedtime", so we modify the Timing Struc to raise the "every_flag" flag.
        In the second case, it means "every night", and we create a new Timing Struc for it.
        We could instead add the "every day" Freq Struc, but prefer to leave those deductions to later reasoning.
        """

        groupname2group = trim_dictionary(match_object.groupdict())

        if 'timing' in groupname2group:
            timing_start = match_object.start('timing')
            timing = parse.position2struc(timing_start)
            if timing.landmark == 'bedtime':
                timing.every_flag = True
                return [timing]
            else:
                # we have no idea what is going on. Fail this rule.
                return None
        else:
            # saying "do something nightly" at the end of the sentence just means "do it every day at night"
            timing = Timing(landmark = 'night', relation = 'at', typ = 'day_part', offset = None, every_flag = True)
            return [timing]

    rule = Rule_ExtractStrucs(   name = 'nightly',
                                    search_patterns = [pattern],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule

def rule_assembly_indication_pain():

    pattern0 = re.compile(r'''      # "pain in" precedes the target (e.g. "joint"). Both target is mandatory.
                            (?<![a-z])
                            (?P<for>
                                (
                                    to \s (?P<verb>reduce|alleviate)    |
                                    for
                                )
                                (\s the)?
                                \s
                            )?
                            (mild \s to \s moderate \s)?
                            (mild \s)?
                            (moderate \s)?
                            (chronic \s)?
                            (?P<severity>
                                (moderate \s? (to|\-) \s? severe \s)       |
                                severe? \s
                            )?
                            (?P<pain_before_target>                 # mandatory
                                pain e? \s in (\s the)? \s
                            )
                            (
                                (?P<pain_type>                      # mandatory
                                    (?P<stomach>
                                        abdom(i|e)n(al)?    |
                                        abnominal           |
                                        stomach
                                    )                           |
                                    (?P<back>
                                        (low(er)? \s)? back (\s (or|and|AND_CONJ) \s joint)?   |
                                        lumbar (\s disc)?
                                    )                           |
                                    (?P<arthritis>
                                        (osteo)? arthritis
                                    )                           |
                                    (?P<bladder>
                                        bladder         |
                                        dysuria         |
                                        urine (\.|$)    |
                                        urology (\.|$)  |
                                    )?                          |
                                    (?P<breakthrough>
                                        break (\s|\-)? through |
                                        breakthru
                                    )                           |
                                    (?P<gout>gout)              |
                                    (?P<joint>joints?)          |
                                    (?P<knee>knees?)            |
                                    (?P<leg>legs?)              |
                                    (?P<chest>chest)            |
                                    (?P<ear>ears?)              |
                                    (?P<nerve>nerves?)
                                )
                            )
                             (?![a-z])
                         ''', re.X)

    pattern1 = re.compile(r'''
                            (?<![a-z])
                            (?P<for>                                           # mandatory "for ... pain" or "to reduce/alleviate ... pain" unless this is the very end of string, e.g. "apply 2 daily pain."
                                (
                                    to \s (?P<verb>reduce|alleviate)    |
                                    for
                                )
                                \s
                            )?
                            (mild \s to \s moderate \s)?
                            (mild \s)?
                            (moderate \s)?
                            (chronic \s)?
                            (?P<severity>
                                (moderate \s to \s severe \s)       |
                                severe? \s
                            )?
                            (
                                (?P<pain_type>
                                    (?P<stomach>
                                        abdom(i|e)n(al)?    |
                                        stomach
                                    )                           |
                                    (?P<back>
                                        (low(er)? \s)? back (\s (or|and|AND_CONJ) \s joint)?   |
                                        lumbar (\s disc)?
                                    )                           |
                                    (?P<arthritis>
                                        (osteo)? arthritis
                                    )                           |
                                    (?P<bladder>bladder)        |
                                    (?P<breakthrough>
                                        break (\s|\-)? through |
                                        breakthru
                                    )                           |
                                    (?P<gout>gout)              |
                                    (?P<joint>joints?)          |
                                    (?P<knee>knees?)            |
                                    (?P<leg>legs?)              |
                                    (?P<chest>
                                        chest       |
                                        heart
                                    )                           |
                                    (?P<ear>ears?)              |
                                    (?P<nerve>nerves?)
                                ) \s?
                            )?
                            (?P<pain_after_target>                  # mandatory
                                 paine      |
                                 pains      |
                                 pain       |
                                 paion      |
                                 paon       |
                                 pian       |
                                 aches?
                            )
                             (?! \s or \s fever)
                             (?![a-z])
                         ''', re.X)



    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        """
        for pain indication (e.g. "for shoulder pain", "to alleviate severe chest pain"
        """

        groupname2group = trim_dictionary(match_object.groupdict())

        # either pain_before_target or pain_after_target is mandatory


        # "for" is mandatory unless at end of sentence or is part of a list of indications ("for heartburn and pain"). E.g. we allow "as needed pain.", but only at the end.
        if 'for' not in groupname2group and 'pain_before_target' not in groupname2group:
            right_context = trim(right_context)
            not_end_of_sentence = right_context and right_context != ' ' and right_context[0] not in ('.', ',', '*')
            left_context = trim(left_context)
            and_conj_immediately_to_the_left = re.compile(r'AND_CONJ\s?$')
            is_part_of_indication_sequence = and_conj_immediately_to_the_left.search(left_context)  and 'INDICATION' in left_context

            if not_end_of_sentence and not is_part_of_indication_sequence:
                return None

            to_is_immediately_to_the_left = re.compile(r'to\s?$')      # Even at end of sentence, avoid labeling "apply to rash" as Indication, not Site.
            if to_is_immediately_to_the_left.search(left_context):
                return None

        if 'severity' in groupname2group:
            severity = True
        else:
            severity = False

        if 'verb' in groupname2group:
            verb = groupname2group['verb']
            if verb in ('reduce', 'alleviate'):
                verb = 'reduce'
        else:
            verb = None

        if 'back' in groupname2group:
            pain_type = 'back'
        elif 'breakthrough' in groupname2group:
            pain_type = 'breakthrough'
        elif 'chest' in groupname2group:
            pain_type = 'chest'
        elif 'ear' in groupname2group:
            pain_type = 'ear'
        elif 'joint' in groupname2group:
            pain_type = 'joint'
        elif 'knee' in groupname2group:
            pain_type = 'knee'
        elif 'nerve' in groupname2group:
            pain_type = 'nerve'
        elif 'stomach' in groupname2group:
            pain_type = 'stomach'
        elif 'arthritis' in groupname2group:        # arthritis pain is just arthritis, not a type of pain.
            indication = Indication(condition = 'arthritis', constituents = match_object.group())
            return [indication]
        elif 'bladder' in groupname2group:        # bladder pain should not be transduced as "take for bladder. Also: take for pain".
            indication = Indication(condition = 'bladder', constituents = match_object.group())
            return [indication]
        elif 'gout' in groupname2group:        # see "bladder pain" above
            indication = Indication(condition = 'gout', constituents = match_object.group())
            return [indication]
        elif 'leg' in groupname2group:        # see "bladder pain" above
            indication = Indication(condition = 'leg', constituents = match_object.group())
            return [indication]
        elif 'pain_type' in groupname2group:
            pain_type = groupname2group['pain_type']
        else:
            pain_type = ''



        struc = IndicationPain(pain_type, severity, verb, constituents = match_object.group())
        return [struc]

    rule = Rule_ExtractStrucs(   name = 'pain indication',
                                    search_patterns = [pattern0, pattern1],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule

def rule_assembly_indication():

    pattern0 = re.compile(r'''              # There are 3 patterns only because of the Python limitation of 100 named groups
                                (?<![a-z])
                                (?P<for>for \s?)?               # "for" is mandatory unless at end of sentence or preceded by "INDICATION AND "
                                                                # (e.g. "for heart and blood pressure"), as enforced by
                                                                # replacement_proc. E.g. we allow "as needed pain.", but only at the end.
                                (the \s)?
                                (severe? \s)?
                                (chronic \s)?
                                (?P<indication>
                                    (?P<acid_reflux>
                                        (acid \s)? reflux       |
                                        acid \s (in \s)? stomach|
                                        acid \s?\.?$            |
                                        gerd
                                    )                               |
                                    (?P<acne>
                                        acne
                                    )                               |
                                    (?P<agitation>
                                        agg?itation
                                    )                               |
                                    (?P<allergies>
                                        (severe \s)?
                                        (
                                            all?erg(y|ies)          |
                                            allergic \s symptoms?   |
                                            allergic \s reactions?
                                        )
                                    )                               |
                                    (?P<anemia>anemia
                                    )                               |
                                    (?P<anxiety>
                                        anxiety                 |
                                        anxious \s feelings?
                                    )                               |
                                    (?P<arthritis>
                                        arthritis (\s pain)?
                                    )                               |
                                    (?P<asthma>
                                        asthma  |
                                        sthma
                                    )                               |
                                    (?P<bladder>
                                        bladder
                                        (\s (spasms?|pain))?
                                    )                               |
                                    (?P<blood_pressure>
                                        (
                                            ((to \s)? control \s)?
                                            (your \s)?
                                            (systolic \s)?
                                            (?<!low\s)(?<!high\s)
                                            blood \s pressure
                                        )                                       |
                                        blood \s (perssure|preassure|presure)   |
                                        bloodpressure                           |
                                        systolic \s pressure                    |
                                        pressure \s over \s QUANT               |  # e.g. "for pressure over 160"
                                        bp                                      |
                                        blood \s? (?=\.|$)                                  # "take one tablet daily for blood" occurs very frequently, but I am not sure this really refers to blood pressure.
                                    )                               |
                                    (?P<high_blood_pressure>
                                        hypertension                            |
                                        hypetension                             |
                                        hypertencion                            |
                                        (
                                            ((to \s)? control \s)?
                                            (your \s)?
                                            (high \s)
                                            (systolic \s)?
                                            (?<!low\s)
                                            blood \s pressure
                                        )                                       |
                                        pressure \s over \s QUANT                   # e.g. "for pressure over 160"
                                    )                               |
                                    (?P<blood_circulation>
                                        (blood \s)? circulation
                                    )                               |
                                    (?P<blood_sugar>
                                        (blood \s)? sugar
                                    )                               |
                                    (?P<blood_thinner>
                                        blood \s thinners?  |
                                        bloodthinner        |
                                        thinner             |
                                        anticoagulants?
                                    )                               |
                                    (?P<bones>
                                        bone s?
                                    )                               |
                                    (?P<breathing>
                                        breathing (\s problems?)?           |
                                        (?<!of \s) breathe? \s? (?=\.|$)
                                    )                               |
                                    (?P<calcium>
                                        calcium
                                        (\s supplements?)?
                                    )                               |
                                    (?P<chest_pain>
                                        chest \s discomfort |
                                        angina
                                    )                               |
                                    (?P<cholesterol>
                                            (high \s)?
                                            (
                                                cholesterol     |
                                                cholsterol      |
                                                cholestserol    |
                                                cholestrol      |
                                                cholestriol     |
                                                cholestorol     |
                                                cholestoral     |
                                                cholestesrol    |
                                                cholesteroly    |
                                                cholesterole    |
                                                cholestero      |
                                                choesterol      |
                                                choelsterol     |
                                                choelesterol    |
                                                chloestrol
                                            )
                                    )                               |
                                    (?P<concentration>
                                        concentration
                                    )                               |
                                    (?P<congestion>
                                        congestion
                                    )                               |
                                    (?P<constipation>
                                        constipation        |
                                        costipation         |
                                        cosntipation        |
                                        constpation         |
                                        constopation        |
                                        constiption         |
                                        constipationaily    |
                                        constipatio         |
                                        constipatino        |
                                        constipatin         |
                                        constipati          |
                                        constiopation       |
                                        constiaption        |
                                        constapation        |
                                        stool \s softener   |
                                        stool \s softner    |
                                        (?<!after\seach\s) bowel \s movement
                                    )                               |
                                    (?P<cough_wheeze>
                                        cough(ing)? \/ wheez(e|ing)                     |
                                        cough(ing)? \s (and|AND_CONJ) \s wheez(e|ing)   |
                                        wheezing \s (and|AND_CONJ) \s cough(ing)?
                                    )                               |
                                    (?P<cough_congestion>
                                        cough (ing)? \s (and|AND_CONJ|or) \s congestion
                                    )                               |
                                    (?P<cough>
                                        cough (ing)?
                                    )                               |
                                    (?P<cramps>
                                        (?<!muscle\s)
                                        (?<!stomach\s)
                                        (?<!abdominal\s)
                                        (?<!abdomen\s)
                                        (menstrual \s)?
                                        cramp(s|ing)?
                                    )                               |
                                    (?P<dementia>
                                        dementia
                                    )                               |
                                    (?P<depression>
                                        depression      |
                                        depresion
                                    )                               |
                                    (?P<diabetes>
                                        diabetes (\s mell?itus)?|
                                        diabeties               |
                                        diabetis                |
                                        diabets
                                    )                               |
                                    (?P<dialysis>
                                        dialysis
                                    )                               |
                                    (?P<diarrhea>
                                        diarrhea    |
                                        diarrehea   |
                                        diarreha    |
                                        diarrhae    |
                                        diarhea     |
                                        diarrea     |
                                        diarrhoea
                                    )                               |
                                    (?P<digestion>
                                        digestion
                                    )                               |
                                    (?P<dizziness>
                                        dizziness?      |
                                        diziness?       |
                                        dizz?y(ness)?
                                    )                               |
                                    (?P<dry_eyes>
                                        dry \s (itchy \s)? eyes?
                                    )                               |
                                    (?P<dry_skin>
                                        dry (\s itchy)? \s skin
                                    )                               |
                                    (?P<dryness>
                                        dryness
                                    )
                                )
                                (\s symptoms?)?
                                (\s discomfort?)?
                                (?![a-z])
                            ''', re.X)

    pattern1 = re.compile(r'''
                                (?<![a-z])
                                (?P<for>for \s?)?               # "for" is mandatory unless at end of sentence or preceded by "INDICATION AND "
                                                                # (e.g. "for heart and blood pressure"), as enforced by
                                                                # replacement_proc. E.g. we allow "as needed pain.", but only at the end.
                                (the \s)?
                                (severe? \s)?
                                (chronic \s)?
                                (?P<indication>
                                    (?P<eczema>
                                        eczema
                                    )                               |
                                    (?P<edema>                              # edema is the same as "swelling", but it's a different register, so we use both.
                                        (edema|adema)
                                        (\s?\-? \(? swelling \)? )?
                                    )                               |
                                    (?P<fever>
                                        fevere?
                                    )                               |
                                    (?P<fungal_infection>
                                        fungal \s infections?   |
                                        fungus (\s infections?)?
                                    )                               |
                                    (?P<gas>
                                        gas         |
                                        flatulence  |
                                        bloat(ing)?
                                    )                               |
                                    (?P<glaucoma>glaucoma)          |
                                    (?P<gout>gout)                  |
                                    (?P<headache>head \s? aches?)   |
                                    (?P<heartburn>
                                        heart (\-|\s)? burn
                                    )                               |
                                    (?P<heartrate>
                                        ((to \s)? control \s)?
                                        (
                                            heart (\-|\s)? rate     |
                                            heart \s? beat          |
                                            heart \s (and|AND_CONJ) \s palpitations?    |
                                            palpitations?
                                        )
                                    )                               |
                                    (?P<heart>
                                        (congestive \s)?
                                        (
                                            heart       |
                                            heaert      |
                                            hearat      |
                                            hearet
                                        )
                                        (\s disease | \s failure | \s pain)?
                                    )                               |
                                    (?P<hemorrhoids>
                                        hemorrhoids?    |
                                        hemm?orr?h?oids?
                                    )                               |
                                    (?P<hiccups>
                                        hiccups?
                                    )                               |
                                    (?P<incontinence>
                                        (urinary \s) ? incontinence     |
                                        urinary \s control      |
                                        urine \s controll?      |
                                        urination \s controll?  |
                                        urinary \s frequency
                                    )                               |
                                    (?P<indigestion>
                                        indigestion     |
                                        dyspepcia       |
                                        dyspepsia
                                    )                               |
                                    (?P<infection>
                                        (?<!tract \s)                   # avoid "urinary tract infection" which has it's own entry
                                        (?<!fungal \s)
                                        infection
                                    )                               |
                                    (?P<inflammation>
                                        inflamm?ation   |
                                        inflamm?acion
                                    )                               |
                                    (?P<insomnia>
                                        sleep\/insomnia     |
                                        to \s sleep         |
                                        sleep \s interr?uptions?    |
                                        sleeping            |
                                        sleep               |
                                        insomnia
                                    )                               |
                                    (?P<iron>
                                        iron
                                        (\s suppl(e|i)ments?)?
                                    )                               |
                                    (?P<irritation>
                                        (to \s reduce \s)?
                                        irritations?
                                    )                               |
                                    (?P<itchy_eyes>
                                        itch(y|ing) \s eyes?
                                    )                               |
                                    (?P<itchy_skin>
                                        itch(y|ing) \s skin
                                    )                               |
                                    (?P<itching>
                                        (
                                            itching     |
                                            itchign     |
                                            itchiness   |
                                            itch y?     |
                                            pruritis    |
                                            pruritus
                                        )
                                        (?!\s (eye|skin))
                                    )
                                )
                                (\s symptoms?)?
                                (\s discomfort?)?
                                (?![a-z])
                            ''', re.X)

    pattern2 = re.compile(r'''
                                (?<![a-z])
                                (?P<for>for \s?)?               # "for" is mandatory unless at end of sentence or preceded by "INDICATION AND "
                                                                # (e.g. "for heart and blood pressure"), as enforced by
                                                                # replacement_proc. E.g. we allow "as needed pain.", but only at the end.
                                (the \s)?
                                (severe? \s)?
                                (chronic \s)?
                                (?P<indication>
                                    (?P<leg_swelling>legs? \s swelling|swelling \s of \s (the \s)? legs?)  |
                                    (?P<leg>legs?)                          |
                                    (?P<lungs>
                                        lungs?
                                    )                                       |
                                    (?P<memory>
                                        memory (\s loss)?   |
                                        memroy
                                    )                                       |
                                    (?P<migraine>
                                        migraine?s?     |
                                        migranes?
                                    )                                       |
                                    (?P<mood>mood (\s stabilization)? )     |
                                    (?P<muscle_spasms>
                                        muscles? \s?
                                        (
                                            spasms?     |
                                            spams?      |
                                            spasams?    |
                                            cramping    |
                                            cramps?     |
                                            (\.|$)
                                        )
                                    )                                       |
                                    (?P<nausea_and_vomiting>
                                        nausea \s? (and|AND_CONJ|or|,)? \s? vomiting     |
                                        nausea \/ vomiting
                                    )                                       |
                                    (?P<nausea>
                                        nausea (asea)?  |
                                        nasuea  |
                                        nauea   |
                                        nause   |
                                        nzusea
                                    )                                       |
                                    (?P<nervousness>
                                        nervousness|nerve s?
                                    )                                       |
                                    (?P<osteoporosis>
                                        osteoporosis (ss?)?     |
                                        osteopirosis

                                    )                                       |
                                    (?P<panic_attack>
                                        (acute \s)?
                                        (anxiety \s (and|AND_CONJ|or) \s)?
                                        panic
                                        (\s? attacks?)?
                                        (\s disorders?)?
                                    )                                       |
                                    (?P<pain_fever>
                                        pain \s  (and|AND_CONJ|or) \s fevere?      |
                                        fevere? \s (and|AND_CONJ|or) \s pain
                                    )                                       |
                                    (?P<potassium>
                                        potassium
                                        (\s supplements?)?
                                    )                                       |
                                    (?P<prostate>prostate)                  |
                                    (?P<rash>
                                        (skin \s)?
                                        rash (es)?
                                    )                                       |
                                    (?P<rhinitis>
                                        rhinitis
                                    )                                       |
                                    (?P<runny_nose>
                                        (runny|running) \s nose     |
                                        rhinorrhea                  |
                                        rhinorrhoea                 |
                                        rhniorrhea                  |
                                        rhinorhea
                                    )                                       |
                                    (?P<seizures>
                                        seizures?   |
                                        epilepsy
                                    )                                       |
                                    (?P<shortness_of_breath>
                                        shor?t(ness)? \s of \s breath   |
                                        dyspnea                         |
                                        sob
                                    )                                       |
                                    (?P<sore_throat>
                                        sore \s throat
                                    )                                       |
                                    (?P<spasms>
                                        spasms?     |
                                        spasams?
                                    )                                       |
                                    (?P<stomach_cramp>
                                        (
                                            stomach         |
                                            abdomin(al)?    |
                                            abdomen(al)?
                                        )
                                        \s
                                        (cramp(ing)?|spasm)s?
                                    )                                       |
                                    (?P<stomach>
                                        (
                                            stomach e?      |
                                            stamach         |
                                            stoamch         |
                                            stomahc         |
                                            stoma c?        |
                                            abdomin(al)?    |
                                            abdomen(al)?
                                        )
                                        (?!=\s (cramp|spasm))
                                    )                                       |
                                    (?P<stuffy_nose>
                                        stuffy \s nose          |
                                        nasal \s congestion     |
                                        nasal \s drainage
                                    )                                       |
                                    (?P<swelling>
                                        (?P<reduce_swelling>to \s reduce \s)?    #needed, in part, to match the dictionary
                                        (ankle \s)?
                                        swelling
                                    )                                       |
                                    (?P<thyroid>
                                        thyroids?       |
                                        hypothyroidism  |
                                        hyperthyroidism |
                                        tyroids?
                                    )                                       |
                                    (?P<tremors>
                                        tremors?
                                    )                                       |
                                    (?P<triglycerides>
                                        triglycerides?  |
                                        triglycirides?  |
                                        triglycrides?
                                    )                                       |
                                    (?P<urinary_tract_infection>
                                        urinary \s tract? \s infections? |
                                        urinary \s problems?    |
                                        urinary \s symptoms?    |
                                        urinary \s pain         |
                                        uti
                                    )                                       |
                                    (?P<vertigo>
                                        vertigo
                                    )                                       |
                                    (?P<vomiting>
                                        vomitt?(ing)?
                                    )                                       |
                                    (?P<wheezing_cough>
                                        wheezing \s cough
                                    )                                       |
                                    (?P<wheezing>
                                        wheezing    |
                                        wheeze      |
                                        wheez       |
                                        whezcgh     |
                                        whezzing
                                    )
                                )
                                (\s symptoms?)?
                                (\s discomfort?)?
                                (?![a-z])
                            ''', re.X)


    pattern3 = re.compile(r'''
                                                              # "pain or fever", "wheezing cough" may be 2 indications already processed
                                (?P<sequence_of_indications>
                                    (?P<indication_0>INDICATION)
                                    \s
                                    (?P<separate_indications>
                                        (and|AND_CONJ|or) \s
                                    )?
                                    (?P<indication_1>INDICATION)
                                )
                            ''', re.X)



    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        """
        for indication (e.g. "for wheezing cough")
        """

        groupname2group = trim_dictionary(match_object.groupdict())

        # "for" is mandatory unless at end of sentence or is part of a list of indications ("for pain and heartburn"). E.g. we allow "as needed pain.", but only at the end.
        if 'for' not in groupname2group and 'reduce_swelling' not in groupname2group and 'sequence_of_indications' not in groupname2group:
            right_context = trim(right_context)
            not_end_of_sentence = right_context and right_context != ' ' and right_context[0] not in ('.', ',', '*')
            left_context = trim(left_context)
            and_conj_immediately_to_the_left = re.compile(r'AND_CONJ\s?$')
            is_part_of_indication_sequence = and_conj_immediately_to_the_left.search(left_context)  and 'INDICATION' in left_context

            if not_end_of_sentence and not is_part_of_indication_sequence:
                return None

            to_is_immediately_to_the_left = re.compile(r'to\s?$')      # Even at end of sentence, avoid labeling "apply to rash" as Indication, not Site.
            if to_is_immediately_to_the_left.search(left_context):
                return None

        if 'acid_reflux' in groupname2group:
            indication = 'acid reflux'
        elif 'acne' in groupname2group:
            indication = 'acne'
        elif 'agitation' in groupname2group:
            indication = 'agitation'
        elif 'allergies' in groupname2group:
            indication = 'allergies'
        elif 'anemia' in groupname2group:
            indication = 'anemia'
        elif 'anxiety' in groupname2group:
            indication = 'anxiety'
        elif 'arthritis' in groupname2group:
            indication = 'arthritis'
        elif 'asthma' in groupname2group:
            indication = 'asthma'
        elif 'bladder' in groupname2group:
            indication = 'bladder'
        elif 'high_blood_pressure' in groupname2group:
            indication = 'high blood pressure'
        elif 'blood_pressure' in groupname2group:
            indication = 'blood pressure'
        elif 'blood_circulation' in groupname2group:
            indication = 'blood circulation'
        elif 'blood_sugar' in groupname2group:
            indication = 'diabetes'
        elif 'blood_thinner' in groupname2group:
            indication = 'blood thinner'
        elif 'bones' in groupname2group:
            indication = 'bones'
        elif 'breathing' in groupname2group:
            indication = 'breathing'
        elif 'calcium' in groupname2group:
            indication = 'calcium'
        elif 'chest_pain' in groupname2group:
            indication = IndicationPain(pain_type = 'chest', severity = False, verb = None, constituents = match_object.group())
            return [indication]
        elif 'cholesterol' in groupname2group:
            indication = 'cholesterol'
        elif 'concentration' in groupname2group:
            indication = 'concentration'
        elif 'congestion' in groupname2group:
            indication = 'congestion'
        elif 'constipation' in groupname2group:
            indication = 'constipation'
        elif 'cough_congestion' in groupname2group:
            indication = 'cough and congestion'
        elif 'cough' in groupname2group:
            indication = 'cough'
        elif 'cough_wheeze' in groupname2group:
            indication = 'cough/wheeze'
        elif 'cramps' in groupname2group:
            indication = 'cramps'
        elif 'dementia' in groupname2group:
            indication = 'dementia'
        elif 'depression' in groupname2group:
            indication = 'depression'
        elif 'diabetes' in groupname2group:
            indication = 'diabetes'
        elif 'digestion' in groupname2group:
            indication = 'digestion'
        elif 'dialysis' in groupname2group:
            indication = 'dialysis'
        elif 'diarrhea' in groupname2group:
            indication = 'diarrhea'
        elif 'dizziness' in groupname2group:
            indication = 'dizziness'
        elif 'dryness' in groupname2group:
            indication = 'dryness'
            if 'eye' in left_context:
                indication = 'dry eyes'
            elif 'skin' in left_context:
                indication = 'dry skin'
            else:
                left_context_pattern = re.compile('(?P<site>SITE)')
                site_obj = left_context_pattern.search(left_context)
                if site_obj:
                    site_start = site_obj.start('site')
                    site = parse.position2struc(site_start)
                    if site.value == 'eye':
                        indication = 'dry eyes'
                    elif site.value == 'skin':
                        indication = 'dry skin'
        elif 'dry_eyes' in groupname2group or 'itchy_eyes' in groupname2group:
            indication = 'dry eyes'
        elif 'dry_skin' in groupname2group:
            indication = 'dry skin'
        elif 'eczema' in groupname2group:
            indication = 'eczema'
        elif 'edema' in groupname2group:
            indication = 'edema'
        elif 'fever' in groupname2group:
            indication = 'fever'
        elif 'fungal_infection' in groupname2group:
            indication = 'fungal infection'
        elif 'gas' in groupname2group:
            indication = 'gas'
        elif 'glaucoma' in groupname2group:
            indication = 'glaucoma'
        elif 'gout' in groupname2group:
            indication = 'gout'
        elif 'headache' in groupname2group:
            indication = 'headache'
        elif 'heartburn' in groupname2group:
            indication = 'heartburn'
        elif 'heartrate' in groupname2group:
            indication = 'heart rate'
        elif 'heart' in groupname2group:
            indication = 'heart'
        elif 'hemorrhoids' in groupname2group:
            indication = 'hemorrhoids'
        elif 'hiccups' in groupname2group:
            indication = 'hiccups'
        elif 'incontinence' in groupname2group:
            indication = 'incontinence'
        elif 'indigestion' in groupname2group:
            indication = 'indigestion'
        elif 'infection' in groupname2group:
            indication = 'infection'
        elif 'inflammation' in groupname2group:
            indication = 'inflammation'
        elif 'insomnia' in groupname2group:
            indication = 'sleep/insomnia'
        elif 'iron' in groupname2group:
            indication = 'iron'
        elif 'irritation' in groupname2group:
            indication = 'irritation'
        elif 'itching' in groupname2group:
            indication = 'itching'
        elif 'itchy_eyes' in groupname2group:
            indication = 'dry eyes'
        elif 'itchy_skin' in groupname2group:
            indication = 'dry skin'
        elif 'leg_swelling' in groupname2group:
            indication = 'leg swelling'
        elif 'leg' in groupname2group:
            indication = 'leg'
        elif 'lungs' in groupname2group:
            indication = 'lungs'
        elif 'memory' in groupname2group:
            indication = 'memory'
        elif 'migraine' in groupname2group:
            indication = 'migraine'
        elif 'mood' in groupname2group:
            indication = 'mood'
        elif 'muscle_spasms' in groupname2group:
            indication = 'muscle spasms'
        elif 'nausea' in groupname2group:
            indication = 'nausea'
        elif 'nausea_and_vomiting' in groupname2group:
            indication = 'nausea and vomiting'
        elif 'nervousness' in groupname2group:
            indication = 'nervousness'
        elif 'osteoporosis' in groupname2group:
            indication = 'osteoporosis'
        elif 'pain_fever' in groupname2group:
            indication = 'pain or fever'
        elif 'panic_attack' in groupname2group:
            indication = 'panic attack'
        elif 'potassium' in groupname2group:
            indication = 'potassium'
        elif 'prostate' in groupname2group:
            indication = 'prostate'
        elif 'rash' in groupname2group:
            indication = 'rash'
        elif 'rhinitis' in groupname2group:
            indication = 'rhinitis'
        elif 'runny_nose' in groupname2group:
            indication = 'runny nose'
        elif 'seizures' in groupname2group:
            indication = 'seizures'
        elif 'shortness_of_breath' in groupname2group:
            indication = 'shortness of breath'
        elif 'sore_throat' in groupname2group:
            indication = 'sore throat'
        elif 'spasms' in groupname2group:
            indication = 'spasms'
        elif 'stomach_cramp' in groupname2group:
            indication = 'stomach cramps'
        elif 'stomach' in groupname2group:
            indication = 'stomach'
        elif 'stuffy_nose' in groupname2group:
            indication = 'stuffy nose'
        elif 'swelling' in groupname2group:
            indication = 'swelling'
        elif 'thyroid' in groupname2group:
            indication = 'thyroid'
        elif 'tremors' in groupname2group:
            indication = 'tremors'
        elif 'triglycerides' in groupname2group:
            indication = 'triglycerides'
        elif 'urinary_tract_infection' in groupname2group:
            indication = 'urinary tract infection'
        elif 'vertigo' in groupname2group:
            indication = 'vertigo'
        elif 'vomiting' in groupname2group:
            indication = 'vomiting'
        elif 'wheezing_cough' in groupname2group:
            indication = 'wheesing cough'
        elif 'wheezing' in groupname2group:
            indication = 'wheezing'

        elif 'sequence_of_indications' in groupname2group:
            indication_0_start = match_object.start('indication_0')
            indication_0 = parse.position2struc(indication_0_start)

            indication_1_start = match_object.start('indication_1')
            indication_1 = parse.position2struc(indication_1_start)

            separate_indications = True if 'separate_indications' in groupname2group else False

            if set([indication_0.condition, indication_1.condition]) == set(['pain', 'fever']):
                indication = 'pain or fever'
            elif set([indication_0.condition, indication_1.condition]) == set(['wheezing', 'cough']) and not separate_indications:
                # "wheezing and cough" != "wheezing cough", hence we conditiong this on separate_indications
                indication = 'wheesing cough'
            elif set([indication_0.condition, indication_1.condition]) == set(['nausea', 'vomiting']):
                indication = 'nausea and vomiting'
            elif set([indication_0.condition, indication_1.condition]) == set(['stomach', 'cramps']) and not separate_indications:
                # "pain in the stomach and cramps" != "stomach cramps", hence we conditiong this on separate_indications
                indication = 'stomach cramps'
            else:
                return [indication_0, Struc(', ', []), indication_1]
        indication = trim(indication)

        struc = Indication(condition = indication, constituents = match_object.group())

        return [struc]


    rule = Rule_ExtractStrucs(   name = 'indication',
                                    search_patterns = [pattern0, pattern1, pattern2, pattern3],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule


def rule_assembly_apply_to_site():

    pattern = re.compile(r'''
                                (?<![a-z])
                                (                                           # mandatory for most sites
                                    (?P<relation>to|onto|on)
                                    \s
                                )?
                                (the \s)?
                                (?<!for \s)
                                (?P<site>
                                    (?P<acne_areas>
                                        acne
                                        (\s areas?)?
                                        (\s on \s face)?
                                    )                                               |
                                    (?P<back>
                                        (?<!lower\s)
                                        (?<!low\s)
                                        back
                                    )                                               |
                                    (?P<body>
                                        (whole \s)?
                                        (entire \s)?
                                        (upper \s)?
                                        body
                                    )                                               |
                                    (?P<cheek>
                                        cheeks?
                                    )                                               |
                                    (?P<diaper_area>
                                        diapers? \s areas?
                                    )                                               |
                                    (?P<eyelash>
                                        (
                                            upper \s (and|AND_CONJ) \s lower    |
                                            both
                                        )?
                                        \s?
                                        (
                                            eyelash(es)?    |
                                            lashes          |
                                            eye \s lash(es)?
                                        )
                                        (\s line)?          #  "apply to lashes line"
                                        (?P<eyelashes_both_eyes>
                                            \s
                                            ((of|in|on) \s)?
                                            (both|each) \s
                                            eyes?
                                        )?
                                    )                                               |
                                    (?P<face>
                                        face
                                    )                                               |
                                    (?P<feet>
                                        ( (both|bilateral) \s)?
                                        (affected \s)?
                                        (
                                            feet
                                        )
                                    )                                               |
                                    (?P<foot>
                                        (affected \s)?
                                        (right \s)?
                                        (left \s)?
                                        foot
                                    )                                               |
                                    (?P<hair>
                                        hair
                                    )                                               |
                                    (?P<hands>
                                        ((both|bilateral) \s)?
                                        hands?
                                    )                                               |
                                    (?P<joint>
                                        ((affected|painful)? \s)?
                                        joints?
                                    )                                               |
                                    (?P<knees>
                                        ( (both|bilateral) \s)?
                                        (affected \s)?
                                        knees?
                                    )                                               |
                                    (?P<legs>
                                        (affected \s)?
                                        legs?
                                    )                                               |
                                    (?P<lesion>
                                        (skin \s)?
                                        lesions?
                                    )                                               |
                                    (?P<lower_back>
                                        lower \s back   |
                                        low \s back
                                    )                                               |
                                    (?P<mouth>
                                        (by \s)?                        # e.g "rinse by mouth" which really means "rinse the mouth"
                                        (in \s)?                        # eg. "rinse in the mouth". But could be route = oral, e.g. "dissolve in mouth"
                                        (the \s)?
                                        mouth
                                    )                                               |
                                    (?P<nails>
                                        (all \s)?
                                        (affected \s)?
                                        (toe \s?)?
                                        nails?
                                    )                                               |
                                    (?P<painful_area>
                                        (most \s)?
                                        (
                                            painful \s areas?   |
                                            pain \s areas?      |
                                            area \s of \s pain  |
                                            affected \s pain \s area s?
                                        )
                                    )                                               |
                                    (?P<rash>
                                        rash (\s areas?)?   |
                                        areas? \s of \s rash
                                    )                                               |
                                    (?P<scalp>
                                        scalp
                                    )                                               |
                                    (?P<skin>
                                        (clean \s)?
                                        ((affected \s)? dry \s)?
                                        (itchy \s)?
                                        skin
                                        (\s areas?)?
                                    )                                               |
                                    (?P<vagina>
                                        (to|per)                            # mandatory
                                        (the \s)?
                                        \s
                                        (
                                            vaginal \s opening  |
                                            vagina
                                        )
                                    )                                               |
                                    (?P<warts>
                                        warts?
                                    )                                               |
                                    (?P<wound>
                                        wounds?
                                    )                                               |
                                    (?P<affected_area>
                                        affected \s area (s|\(s\))?         |
                                        affected \s skin \s area (s|\(s\))? |
                                        affected \s skin                    |
                                        itch(y|ing)? \s area s?             |
                                        infected \s area s?                 |
                                        injured \s area s?                  |
                                        area
                                    )
                                )
                                (?![a-z])
                            ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        """
        Apply to affected areas
        """

        groupname2group = trim_dictionary(match_object.groupdict())

        if 'affected_area' in groupname2group:
            site = 'affected area'
        elif 'acne_areas' in groupname2group:
            site = 'acne areas'
        elif 'back' in groupname2group:
            site = 'back'
        elif 'body' in groupname2group:
            site = 'body'
        elif 'cheek' in groupname2group:
            site = 'cheek'
        elif 'eyelash' in groupname2group:
            site = 'eyelash'
        elif 'face' in groupname2group:
            site = 'face'
        elif 'feet' in groupname2group:
            site = 'feet'
        elif 'foot' in groupname2group:
            site = 'foot'
        elif 'diaper_area' in groupname2group:
            site = 'diaper area'
        elif 'boot' in groupname2group:
            site = 'foot'
        elif 'hair' in groupname2group:
            site = 'hair'
        elif 'hands' in groupname2group:
            site = 'hands'
        elif 'joint' in groupname2group:
            site = 'joint'
        elif 'knees' in groupname2group:
            site = 'knees'
        elif 'legs' in groupname2group:
            site = 'legs'
        elif 'lesion' in groupname2group:
            site = 'lesion'
        elif 'lower_back' in groupname2group:
            site = 'lower back'
        elif 'mouth' in groupname2group:
            site = 'mouth'
        elif 'nails' in groupname2group:
            site = 'nails'
        elif 'painful_area' in groupname2group:
            site = 'painful area'
        elif 'rash' in groupname2group:
            site = 'rash'
        elif 'scalp' in groupname2group:
            site = 'scalp'
        elif 'skin' in groupname2group:
            site = 'skin'
        elif 'vagina' in groupname2group:
            site = 'vagina'
        elif 'warts' in groupname2group:
            site = 'warts'
        elif 'wound' in groupname2group:
            site = 'wound'

        if 'relation' in groupname2group:
            relation = groupname2group['relation']
        elif site in ('affected area', 'diaper area', 'eyelash', 'mouth', 'vagina'):
            # For "affected area" there are very frequent cases that omit relation, e.g. "apply affected area twice a day"
            # "rinse mouth" also has no relation.
            relation = ''
        else:
            return None


        if site == 'mouth':
            # site can be mouth only if directive is "rinse" or "swish"
            # Otherwise it should be route (e.g. "dissolve 1 tablet in mouth"
            if 'rinse' not in left_context and 'swish' not in left_context:
                left_context_pattern = re.compile('(?P<left_dir>DIRECTIVE)')
                left_dir_match_obj = left_context_pattern.search(left_context)
                if not left_dir_match_obj:
                    return None
                left_dir_start = left_dir_match_obj.start('left_dir')
                left_dir = parse.position2struc(left_dir_start)
                if left_dir.value not in ('rinse', 'swish'):
                    return None


        site_struc = Site(value = site, relation = relation, constituents = [match_object.group()])
        return [site_struc]

    rule = Rule_ExtractStrucs(   name = 'apply_to_site',
                                    search_patterns = [pattern],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule

def rule_assembly_instillables_site():

    pattern0 = re.compile(r'''
                                (?<![a-z])
                                (in \s (eye|ear|nostril) \s)?               # there is a bizzare form "instill 1 drop in eye in each eye"
                                ((?P<relation>into|in \s to|in|to) \s)?     # optional (there is common use of "2 sprays each nostril" without preps
                                (the \s)?
                                (?P<which>left|right|affected|surgery|operative|each|both|bilateral|per|(?P<quant>QUANT)) \s         # mandatory
                                (?P<site>
                                    (?P<cheek>
                                        (side \s of \s)? cheek
                                    )                           |
                                    (?P<ear>
                                        ear
                                    )                           |
                                    (?P<eye>
                                        eye(?!lid)
                                    )                           |
                                    (?P<nostril>
                                        nostrill?   |
                                        nostil
                                    )
                                )
                                (s| \( s \) )?
                                (?![a-z])
                            ''', re.X)

    pattern1 = re.compile(r'''                              # no "which". E.g. "2 sprays in nostrils" But is plural.
                                (?<![a-z])
                                ((?P<relation>into|in \s to|in|to) \s)?     # optional (there is common use of "2 sprays each nostril" without preps
                                (the \s)?
                                (?P<site>
                                    (?P<cheek>
                                        cheeks
                                    )                           |
                                    (?P<ear>
                                        ears    |
                                        ear\(s\)
                                    )                           |
                                    (?P<eye>
                                        eyes        |
                                        eye\(s\)
                                    )                           |
                                    (?P<nostril>
                                        nostrill?s  |
                                        nostill?s   |
                                        nostill?s\(s\)
                                    )
                                )
                                (?![a-z])
                            ''', re.X)


    pattern2 = re.compile(r'''              # eyelids, which can be upper/lower/left/right/both and combination of these (both upper eyelids, left upper, etc.)
                                            # We recognize only the dominant ones: left/right if used, then lower/upper, then both
                                (?<![a-z])
                                ((?P<relation>into|in \s to|in|to) \s)?     # optional e.g. "apply both eye lids at bedtime"
                                (the \s)?
                                (?P<which_eyelid>                                   # optional
                                    ((left|right|upper|lower|affected|surgery|operative|each|both) \s?)+
                                )?
                                (?P<eyelid>
                                    (eye)?lid
                                    (?P<plural>s)?
                                )
                                (?![a-z])
                            ''', re.X)


    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        """
        (instill|place|apply|use|inhale) N (drops|sprays) (into|to|in) (left|right|affected|each|both|1) (ear|eye|nostril)
        """

        if 'DOSE' in left_context:                      # A dose should be specified before, e.g. 2 drops in each eye
            pass
        elif re.search('DIRECTIVE\s?$', left_context):          # Or, at least the site should be preceded by the right directive
            directive_start = left_context.rfind('DIRECTIVE')
            directive = parse.position2struc(directive_start)
            dir_value = directive.value
            if dir_value not in ('instill', 'place', 'apply', 'use', 'inhale', 'spray'):
                return None


        groupname2group = trim_dictionary(match_object.groupdict())

        if 'cheek' in groupname2group:
            site = 'cheek'
        elif 'ear' in groupname2group:
            site = 'ear'
        elif 'eyelid' in groupname2group:
            site = 'eyelid'
        elif 'eye' in groupname2group:
            if 'eyelid' in right_context or 'eyelid' in left_context:
                # e.g. "apply in both eyes eye lids at bedtime"
                # But we really need to make sure that there is no SITE on the left or the rigth with with value=eyelid. No time for that.
                return None
            site = 'eye'
        elif 'nostril' in groupname2group:
            site = 'nostril'


        if 'which' in groupname2group:
            which = groupname2group['which']
            if which in ('both', 'bilateral', 'per'):
                # example: "use 2 sprays per nostril"
                which = 'each'
            elif which in ('surgery', 'operative'):
                which = 'affected'
        else:
            which = 'each'

        if 'which_eyelid' in groupname2group:       # pattern2 for eyelids
            which_eyelid = groupname2group['which_eyelid']
            if 'left' in which_eyelid:
                which = 'left'
            elif 'right' in which_eyelid:
                which = 'right'
            elif 'lower' in which_eyelid:
                which = 'lower'
            elif 'upper' in which_eyelid:
                which = 'upper'
            elif 'both' in which_eyelid:
                which = 'each'
            elif 'each' in which_eyelid:
                which = 'each'
            elif 'plural' in groupname2group:
                which = 'each'
            else:
                which = None
        elif 'eyelid' in groupname2group:
            # e.g. "apply eyelids"
            if 'plural' in groupname2group:
                which = 'each'
            else:
                which = None

        if 'quant' in groupname2group:      # used in "2 sprays into one nostril"
            quant_start = match_object.start('quant')
            quant = parse.position2struc(quant_start)
            if quant.value == 1 or quant.num_type == 'var':
                which = 'one'
            else:       # hopefully, quant.value == 2, but we default to "each"
                which = 'each'


        if 'relation' in groupname2group:
            relation = groupname2group['relation']
        else:
            relation = ''

        site_struc = Site(value = site, relation = relation, constituents = [match_object.group()], which = which)
        return [site_struc]

    rule = Rule_ExtractStrucs(   name = 'apply_to_instillables_site',
                                    search_patterns = [pattern0, pattern1, pattern2],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule


def rule_assembly_other_site():

    pattern = re.compile(r'''
                                (?<![a-z])
                                (?P<relation>in) \s (the \s)?
                                (?P<site>
                                   (?P<cheek>cheek)
                                )
                                (?![a-z])
                            ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        """
        Dissolve in the cheek and other misc sites
        """

        left_context_pattern = re.compile('(?P<directive>DIRECTIVE)')
        directive_obj = left_context_pattern.search(left_context)
        if not directive_obj:
            return None
        directive_start = directive_obj.start('directive')
        directive = parse.position2struc(directive_start)

        groupname2group = trim_dictionary(match_object.groupdict())
        site = groupname2group['site']

        if 'cheek' in groupname2group:
            site = 'cheek'
            if directive.value != 'dissolve':
                return None

        relation = groupname2group['relation']

        site_struc = Site(value = site, relation = relation, constituents = [match_object.group()])
        return [site_struc]

    rule = Rule_ExtractStrucs(   name = 'apply_to_site',
                                    search_patterns = [pattern],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule




def rule_assembly_as_directed_1():
    """
        -- Deal with "as directed" type of statements.
    """

    pattern0 = re.compile(r'''             # mandatory: "doctor's instructions" or "as your doctor instructed"
                            (?<![a-z])
                            (?P<or_flag>or \s)?
                            (follow|(?P<verb_override>as \s directed|as \s instructed) \s)?
                            (?P<directive>to \s be \s (used|taken) \s)?
                            (as \s)?(per \s)?
                            (?P<your>your \s)?
                            (?P<authority>doctor|physician|provider|prescriber|m\.?d\.?|dr\.?)(\'?s)?       # mandatory
                            \s
                            (?P<verb>directions?|instructions?|orders?|directed|instructed|prescr?ibed|explained)   # mandatory
                            (?P<exactly>exactly)?
                            (?![a-z])
                            ''', re.X)

    pattern1 = re.compile(r'''             # mandatory: "(as directed|per|follow) packaging (instructions)?"
                            (?<![a-z])
                            (?P<or_flag>or \s)?
                            (?P<exactly>exactly \s)?
                            (?P<verb>as \s directed|as \s instructed|as \s prescr?ibed|as \s explained|as \s stated|(as \s)? per|follow) \s      # mandatory
                            (per \s)?
                            (on \s)?
                            (the \s)?
                            ((attached|enclosed) \s)?
                            (patient\'? \s)?
                            (?P<authority_override>doctor\'? \s)?
                            (?P<authority>package|packaging|pack(et)?|pak|sheet) \s?       # mandatory.
                            (?P<verb_override>(instructions?|directions) \s?)?
                            ((enclosed \s)?(?P<sheet>sheet) \s?)?
                            (attached|enclosed)?
                            (?![a-z])
                            ''', re.X)



    pattern2 = re.compile(r'''             # mandatory: "(as directed|per) (instructions)? on packaging"
                            (?<![a-z])
                            (?P<or_flag>or \s)?
                            (?P<exactly>exactly \s)?
                            (?P<verb_override>as \s directed|as \s instructed|as \s prescr?ibed|as \s explained|as \s stated|(as \s)? per) \s      # mandatory
                            (the \s)?
                            ((attached|enclosed) \s)?
                            (patient\'? \s)?
                            (doctor\'? \s)?
                            ((?P<verb>instructions?|directions?) \s)?
                            (on|inside|in|by) (\s the)? \s                                      # mandatory. "inside" is used for "inside box"
                                                                                                # "by" is used for "as instructed by package"
                            (attached|enclosed \s)?
                            (?P<authority_override>package|packaging|pack(et)?|pak|sheet|box) \s?   # mandatory
                            (attached|enclosed)?
                            (?![a-z])
                            ''', re.X)


    pattern3 = re.compile(r'''             # mandatory: "follow (directions/instructions) from doctor"
                            (?<![a-z])
                            (?P<or_flag>or \s)?
                            (follow \s)                                                                 # mandatory
                            (the \s)?
                            (?P<verb>directions?|instructions?) \s                                      # mandatory
                            (from \s|by \s)?
                            (?P<your>your \s)?
                            (?P<authority>doctor|physician|provider|prescriber|m\.?d\.?|dr\.?)(\'?s)?   # mandatory
                            (?P<exactly>exactly)?
                            (?![a-z])
                            ''', re.X)

    pattern4 = re.compile(r'''             # mandatory: "follow instructions (on packaging)?"
                            (?<![a-z])
                            (?P<or_flag>or \s)?
                            (follow) \s                                             # mandatory
                            (?P<exactly>exactly \s)?
                            (the \s)?
                            ((attached|enclosed) \s)?
                            (patients?(\'s?)? \s)?
                            (?P<authority_override2>(doctors?(\'s?)?|package|packaging|pack(et)?|pak|sheet) \s)?
                            (?P<verb>instructions?|directions?) \s?                # mandatory
                            (
                                (?P<exactly2>exactly \s)?
                                (as \s stated \s)?
                                (on|inside|in|by) (\s the)? \s                          # mandatory if group is used
                                (separate \s)?
                                ((attached|enclosed) \s)?
                                (?P<authority_override>package|packaging|pack(et)?|pak|sheet|box) \s?  # mandatory if group is used
                            )?
                            (attached|enclosed)?
                            (?![a-z])
                            ''', re.X)

    pattern5 = re.compile(r'''             # mandatory: as directed (by doctor)?
                            (?<![a-z])
                            (?P<or_flag>or \s)?
                            (?P<directive>to \s be \s (used|taken) \s)?
                            (?P<exactly>exactly \s)?
                            (as|when) \s                                            # mandatory (e.g. "administer when directed by your physician"
                            (?P<verb>directed|instructed|prescr?ibed|explained)     # mandatory
                            (
                                (\.)?
                                \s
                                (by|per) \s                                         # mandatory if group is used
                                (the \s)?
                                (?P<your>your \s)?
                                (?P<authority>doctor|physician|provider|prescriber|m\.?d\.?|dr\.?)(\'s)?  # mandatory if group is used
                                (\s (instructions?|directions))?
                            )?
                            (?![a-z])
                            ''', re.X)


    patterns = [pattern0, pattern1, pattern2, pattern3, pattern4, pattern5]

    def search_proc(txt, pattern, start = None):
        start_pos = start if start else 0
        match_obj = pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        """ Directive modifiers such as "as directed", "per instructions of md", "follow packaging instructions"

        Extracts, if appropriate, the following attributes of the AsDirected() struc:
        verb            string. The action of the prescribing authority: E.g. directed, instructed, explained, prescribed
        authority       string. E.g. doctor, package, etc.
        exactly_flag    True/False if "exactly" is present, e.g. "use exactly as prescribed"
        your_flag       True/False if "your doctor" is present.
        or_flag         True/False if "or as directed" type of struc is present.
        directive       (Optional) String: The action you are to perform with the prescription (similar to DIRECTIVE): take, use, give, inject, etc.

        The difference between verb, directive, and authority is illustrated here: "take as directed by your doctor":
        verb = 'directed'
        authority = 'doctor'
        directive = 'take'
        """

        groupname2group = trim_dictionary(match_object.groupdict())
        constituents = match_object.group()

        authority = None
        if 'authority' in groupname2group:
            authority = groupname2group['authority']
        if 'authority_override' in groupname2group:
            authority = groupname2group['authority_override']
        if 'authority_override2' in groupname2group:
            authority = groupname2group['authority_override2']

        if authority:
            authority = authority.replace('.', '')
            if authority in ('doctor', 'physician', 'md', 'dr', 'provider', 'prescriber'):
                authority = 'doctor'
            elif 'doctor' in authority:     # e.g. "doctor's"
                authority = 'doctor'
            elif 'pack' in authority or 'pak' in authority or authority == 'box':
                authority = 'packaging'
            elif 'instruc' in authority:
                authority = 'instructions'

        if 'sheet' in groupname2group:
            authority = 'sheet'

        if 'verb' in groupname2group:
            verb = groupname2group['verb']
            #verb = trim(verb.replace('as ', '').replace('per', ''))
            if 'direct' in verb:
                verb = 'directed'
            elif 'instr' in verb:
                verb = 'instructed'
            elif 'presc' in verb:
                verb = 'prescribed'
            elif 'expl' in verb:
                verb = 'explained'
            else:
                verb = None
        else:
            verb = None

        if 'verb_override' in groupname2group:
            verb = groupname2group['verb_override']
            if 'direct' in verb:
                verb = 'directed'
            elif 'instr' in verb:
                verb = 'instructed'
            elif 'presc' in verb:
                verb = 'prescribed'
            elif 'expl' in verb:
                verb = 'explained'
            else:
                verb = None

        if 'your' in groupname2group:
            your_flag = True
        else:
            your_flag = False

        if 'exactly' in groupname2group or 'exactly2' in groupname2group:
            exactly_flag = True
        else:
            exactly_flag = False

        if 'or_flag' in groupname2group:
            or_flag = True
        else:
            or_flag = False

        struc = AsDirected(verb = verb, constituents = [constituents], authority = authority,
                            exactly_flag = exactly_flag, your_flag = your_flag, or_flag = or_flag)

        if 'directive' in groupname2group:
            directive = groupname2group['directive']
            if 'use' in directive:
                directive = 'use'
            elif 'take' in directive:
                directive = 'take'
            elif 'remove' in directive:
                directive = None
            elif directive not in Directive.permissible_values:
                directive = None
            if directive:
                struc.directive = directive

        # Now if there is another Directive on the left, consolidate them
        # because each Schedule can only have 1 AS_Directed.
        # Example: "follow instructions on package with meals and at bedtime as directed"
        left_context_pattern = re.compile('(?P<prev_struc>AS_DIRECTED)')
        found_obj = left_context_pattern.search(left_context)
        if found_obj:
            prev_struc_start = found_obj.start('prev_struc')
            prev_struc = parse.position2struc(prev_struc_start)
            if verb and not prev_struc.verb:
                prev_struc.verb = verb
            if authority and not prev_struc.authority:
                prev_struc.authority = authority
            if exactly_flag and not prev_struc.exactly_flag:
                prev_struc.exactly_flag = exactly_flag
            if your_flag and not prev_struc.your_flag:
                prev_struc.your_flag = your_flag
            if or_flag and not prev_struc.or_flag:
                prev_struc.or_flag = or_flag
            return []


        return [struc]

    rule = Rule_ExtractStrucs(   name = 'as_directed_1',
                                    search_patterns = patterns,
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule


def rule_assembly_as_directed_2():
    """ Handles cases such as "... then take as directed ..." where there is a repeated directive in front of the "as directed" phrase.

    We want to accomodate the following sometimes conflicting requirements:
    -   Treat the directive in "instill as directed" or "take as directed" as a genuine directive because it can be used elsewhere as the
        genuine Directive. E.g. "Take as directed 2 tablets". We don't want to parse it as AS_DIRECTED DOSE without any DIRECTIVE.
    -   Swallow the directive immediately preceding "as_directed" if it there is yet another directive to our left, e.g. "take 2 tabs daily, use as directed"
        Here we want to process it as DIRECTIVE DOSE PERIODICITY AS_DIRECTED and not as DIRECTIVE DOSE PERIODICITY DIRECTIVE AS_DIRECTED
    -   The actual dictionary items should be processed purely as AS_DIRECTED (e.g. "use as directed by your doctor." should be only AS_DIRECTED)
        because otherwise they are represented as "DIRECTIVE AS_DIRECTED" in the dictonary and when matching the sig, "take 2 tabs, use as directed"
        we will never ever literally match "use as directed" to "use as directed" in the dictionary because we will be matching "DIRECTIVE/take AS_DIRECTED"
        to "DIRECTIVE/use AS_DIRECTED". But if the sig actually says "use as directed" it is best to translate it as "use", not as "take" or any other directive
        used in the front of the sig.
    """

    pattern = re.compile(r'(?P<directive>DIRECTIVE) (?P<as_dir>AS_DIRECTED)')

    def search_proc(txt, pattern, start = None):
        start_pos = start if start else 0
        match_obj = pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())

        # Skip this rule if there is no other Directive to our left anywhere or immmediately on our right.
        # The reason for skipping is that otherwise we seem to have 2 directives in the sig, and we interpret the onset of
        # the second directive as the onset of a new Admin_Event. E.g. "use as directed. Apply to ..." is interpreted as
        # Use as Directed. Also Apply ..."
        # However, if the sig is nothing but "DIRECTIVE AS_DIRECTED" (ie. left and rigth context are empty), then proceed with
        # this rule because the sig really just says AS_DIRECTED and the Directive in it has no other function outside of AS_DIRECTED.
        #
        # However, even if we are skipping this rule (i.e. not removing the Directive), if the value of AS_DIRECTED.directive is empty,
        # we should fill it with the value from the DIRECTIVE (e.g then "inject as directed by doctor in office" will actually match
        # the dictionary's "inject as directed by doctor" as opposed to "use as directed by doctor".)


        directive_start = match_object.start('directive')
        directive = parse.position2struc(directive_start)

        if directive.value == 'remove':     # "Remove/off" often functions with duration in special ways, so don't swallow it. E.g. "Insert 1 ring for 3 weeks, Then 1 week off as directed"
            return None

        as_dir_start = match_object.start('as_dir')
        as_dir = parse.position2struc(as_dir_start)

        if not as_dir.directive:
            as_dir.directive = directive.value

        left_context_pattern = re.compile('DIRECTIVE')
        previous_directive = left_context_pattern.search(left_context)
        right_context_pattern = re.compile('\W*DIRECTIVE')
        immediately_following_directive = right_context_pattern.match(right_context)
        if not previous_directive and not immediately_following_directive and (left_context or right_context):
            return None


        as_dir.directive = directive.value
        return [as_dir]

    rule = Rule_ExtractStrucs(   name = 'as_directed_2',
                                    search_patterns = [pattern],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule


def rule_assembly_see_instructions():

    pattern0 = re.compile(r'''
                            (?<![a-z])
                            (see |as \s per) \s                                         # mandatory
                            (the \s)?
                            (?P<attached>attached|attachment|attatched|enclosed)? \s?             # in the body of the replacement proc we ensure that
                                                                                        # at least one of <attached> or <authority> is mandatory.
                            ((package|packaging|pack(et)?|pak|rx) \s)?
                            (
                                ((written|printed) \s)?
                                (?P<authority>
                                    label|note|memo|comment     |
                                    direction                   |
                                    instruction (\s sheet)?     |
                                    sheet|information|enclosure|paper|schedule|form
                                )s?
                                \s?
                            )?
                            (for \s (further \s)? (dosing \s)? (instructions|directions))?
                            (?![a-z])
                        ''', re.X)



    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())
        found = match_object.group()

        if 'attached' in groupname2group:
            attached_flag = True
        else:
            attached_flag = False

        if 'authority' in groupname2group:
            authority = groupname2group['authority']
            if authority in ('direction', 'instruction'):
                authority = 'directions'
            elif authority in ('note', 'memo', 'comment'):
                authority = 'note'
            elif authority in ('sheet', 'instruction sheet', 'information', 'enclosure', 'paper', 'schedule', 'form'):
                authority = 'sheet'
        else:
            authority = ''

        if not authority and not attached_flag:
            return None

        struc = SeeInstructions(authority = authority, constituents = [found], attached_flag = attached_flag)


        return [struc]

    rule = Rule_ExtractStrucs(   name = 'see_instructions',
                                    search_patterns = [pattern0],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule


def rule_assembly_as_needed():

    pattern = re.compile(r'''
                                (?<![a-z])
                                (?P<only>only \s)?
                                (                           # mandatory
                                    as          |
                                    when        |
                                    (?P<if>if)
                                )
                                \s
                                (                           # either "needed" or "necessary" are mandatory
                                    (you \s)?
                                    need
                                    (ed|\s it)?
                                |
                                    necessary
                                )
                                (?![a-z])
                            ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())
        found = match_object.group()
        if 'only' in groupname2group:
            only_flag = True
        else:
            only_flag = False

        if 'if' in groupname2group:
            if_flag = True
        else:
            if_flag = False
        dir_struc = AsNeeded(constituents = [found], only_flag = only_flag)
        dir_struc.if_flag = if_flag


        return [dir_struc]

    rule = Rule_ExtractStrucs(   name = 'as_needed',
                                    search_patterns = [pattern],
                                    search_proc = search_proc,
                                    replacement_proc = replacement_proc,
                                    rule_properties = set(['struc_identification']),
                                    change_parse_in_place = True)
    return rule




def rule_assembly_vehicle():

    pattern0 = re.compile(r'''
                            (?<![a-z])
                            (?P<preposition>via|in|with|by|using|per)         # mandatory
                            \s
                            (the \s)?
                            (per \s)?           # 1 vial via per nebulizer
                            (a \s)?
                            (your \s)?
                            (hand (\s|\-) held \s)?
                            (
                                (?P<nebulizer>
                                    nebuli(z|s)er   |
                                    nebulzier       |
                                    neb             |
                                    svn                 # small-volume nebulizer
                                )                                   |
                                (?P<inhaler>
                                    (albuterol \s)? inhaler     |
                                    albuterol                   |
                                    handihaler                  |
                                    mdi                             # metered dose inhaler
                                )
                            )
                            (?![a-z])
                        ''', re.X)

    pattern1 = re.compile(r'''
                                (?P<nebulize_dir>       # "nebulize" as verb is directive inhale via nebulizer
                                    ^nebulize(?=\s)
                                )
                        ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())
        found = match_object.group()
        if 'nebulizer' in groupname2group:
            vehicle = 'nebulizer'
        elif 'nebulize_dir' in groupname2group:
            vehicle = 'nebulizer'
            vehicle = Vehicle(value = vehicle, preposition = '', constituents = [found])
            directive = Directive('inhale')
            return [directive, Struc(' ', []), vehicle]
        elif 'inhaler' in groupname2group:
            vehicle = 'inhaler'
        preposition = groupname2group['preposition']
        struc = Vehicle(value = vehicle, preposition = preposition, constituents = [found])

        return [struc]

    rule = Rule_ExtractStrucs(  name = 'vehicle',
                                search_patterns = [pattern0, pattern1],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)
    return rule






def rule_assembly_varquant():
    """ Identifies and codes Variables from the dictionary. E.g. <<NUM_3>> or <<TIME_0>> or <<DATE_0>> """

    pattern = re.compile(r'''
                            <<
                            (?P<var_type>
                                num     |
                                date    |
                                time
                            )
                            \_
                            (
                                (?P<quant>QUANT)    |
                                (?P<num>\d+)
                            )
                            >>
                        ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())
        var_type = groupname2group['var_type']
        if 'quant' in groupname2group:
            quant_start = match_object.start('quant')
            quant = parse.position2struc(quant_start)
            var_number = quant.value
        elif 'num' in groupname2group:
            var_number = int(groupname2group['num'])
        struc = VarQuant(var_type = var_type, var_number = var_number)
        struc.constituents = ['<<' + struc.value_repr + '>>']

        return [struc]

    rule = Rule_ExtractStrucs(  name = 'dictionary variable identification',
                                search_patterns = [pattern],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)
    return rule






def rule_assembly_maxdose():

    pattern0 = re.compile(r'''                  # Time interval is specified in the "Max" expression, before Quant
                                                # e.g. MDD, or "max daily dose 8 tabs", etc.
                            (?<![a-z])
                            (                   # (optional) prefix "do not take more than"
                                (do \s not | not \s to) \s
                                (
                                    (?P<directive>DIRECTIVE) \s  more \s than  |
                                    exceed (more \s than)?
                                )
                            \s )?
                            (with (\s a)? \s)? (the \s)?
                            (?P<daily_flag>                         # one of the following mdd expressions is mandatory
                                daily \s max(imum)? (\s dose)?                              |
                                max\.?(imum)? \s (daily|per \s day|a \s day) (\s (dose|dosage))?     |
                                max\.?(imum)? \s? (per|a|\/) \s? TIMEUNIT                   |
                                mdd
                            ) \s?
                            (                   # (optional) suffix "do not take more than"
                                (do \s not | not \s to) \s
                                (
                                    (?P<directive2>DIRECTIVE) \s  more \s than  |
                                    exceed (more \s than)?
                                )
                            \s )?
                            ((is|of|\=|\:|\-) \s?)?
                            (                                       # mandatory: Quant (e.g "mdd=5") or Dose ("mdd 5 tabs/24 hours") or Freq ("not more than twice daily")
                                (?P<dose_quant>QUANT)   |
                                (?P<dose>DOSE)          |
                                (?P<freq>FREQ)
                            ) \s?
                            (                                       # optional: "daily" or "24 hours" (even though it is redundant with "daily" above
                                (a|per|in|\/)?\s?
                                (                                   # mandatory if group is chosen:
                                    (day|daily|TIMEUNIT)    |       # include "TIMEUNIT" if "day" was already processed as TIMEUNIT. But then its value has to be "day"
                                    (?P<time_interval>TIMEINTERVAL)
                                )
                            )?
                            (?![a-z])
                            ''', re.X)

    pattern1 = re.compile(r'''                  # Time interval is specified in the TIMEINTERVAL Struc.
                            (?<![a-z])
                            (                   # Mandatory : one of "Maximum" or  "do not take more than" or "do not exceed"
                                (
                                    (do \s not | not \s to) \s
                                    (
                                        (?P<directive>DIRECTIVE) \s  more \s than  |
                                        exceed (more \s than)?
                                    )
                                    \s
                                    ((the \s)? max\.?(imum)? (\s dose)? \s?)?
                                    )
                                |

                                (with (\s a)? \s)? (the \s)? max(imum)? (\s dose)? \s?
                                |
                                not? \s more \s than \s?
                                |
                                limit \s of \s?
                                |
                                until \s reach(ed)? \s?
                            )
                            ((is|of|\=|\:) \s?)?
                            (                                           #mandatory: Quant (e.g "mdd=5") or Dose ("mdd 5 tabs/24 hours") or Freq ("not more than twice daily")
                                (?P<dose_quant>QUANT (\s doses?)?)  |
                                (?P<dose>DOSE)                      |
                                (?P<freq>FREQ)
                            ) \s?
                            (of \s (this \s medicine|ANAPHORA) \s)?
                            (PERIODICITY \s)?                           # There is sometimes redundant "daily", e.g. "do not take more than 1 tab daily per 24/hours"
                            (                                           # mandatory
                                (a|per|in|\/)?\s?
                                (                                       # mandatory if group is chosen:
                                    (?P<daily>daily)                |
                                    (?P<timeunit>TIMEUNIT)          |
                                    (?P<time_interval>TIMEINTERVAL) |
                                    (?P<periodicity>PERIODICITY)        # e.g. maximum 2 tabs every 24 hours
                                )
                            )
                            (?![a-z])
                            ''', re.X)



    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())
        if 'dose' in groupname2group:
            dose_start = match_object.start('dose')
            dose = parse.position2struc(dose_start)
        elif 'dose_quant' in groupname2group:
            dose_quant_start = match_object.start('dose_quant')
            dose_quant = parse.position2struc(dose_quant_start)
            # Try to get it from the left DOSE if possible
            left_context_pattern = re.compile('(?P<dose_on_left>DOSE)')
            left_dose_obj = left_context_pattern.search(left_context)
            form_name = 'unit'
            if left_dose_obj:
                left_dose_start = left_dose_obj.start('dose_on_left')
                left_dose = parse.position2struc(left_dose_start)
                if left_dose.form:
                    form_name = left_dose.form.value
            form = Form(form_name = form_name, plurality = None, constituents = [])
            dose = Dose(quant = dose_quant, form = form)
        elif 'freq' in groupname2group:
            freq_start = match_object.start('freq')
            freq = parse.position2struc(freq_start)
            dose_quant = freq.quant
            # Try to get it from the left DOSE if possible
            left_context_pattern = re.compile('(?P<dose_on_left>DOSE)')
            left_dose_obj = left_context_pattern.search(left_context)
            form_name = 'unit'
            if left_dose_obj:
                left_dose_start = left_dose_obj.start('dose_on_left')
                left_dose = parse.position2struc(left_dose_start)
                if left_dose.form:
                    form_name = left_dose.form.value
            form = Form(form_name = form_name, plurality = None, constituents = [])
            dose = Dose(quant = dose_quant, form = form)
        elif debug:
            raise Exception("Error in parsing max dose in raw sig: ->%s<-\nDose indicator is missing in found expression ->%s<-" %
                            (parse.sig.raw_sig, match_object.group()))
        else:
            return None

        if 'time_interval' in groupname2group:
            time_interval_start = match_object.start('time_interval')
            time_interval = parse.position2struc(time_interval_start)
        elif 'periodicity' in groupname2group:
            periodicity_start = match_object.start('periodicity')
            periodicity = parse.position2struc(periodicity_start)
            time_interval = TimeInterval(quant = periodicity.quant, time_unit = periodicity.time_unit)
        elif 'daily_flag' in groupname2group or 'daily' in groupname2group:     # per 1 day is implied as TimeInterval.
            quant = Quant(num_type = 'int', value = 1)
            if 'daily_flag' in groupname2group:
                time_unit_constituents = match_object.groupdict()['daily_flag']
            else:
                time_unit_constituents = 'daily'
            time_unit = TimeUnit(value = 'day', constituents = [time_unit_constituents])
            time_interval = TimeInterval(quant = quant, time_unit = time_unit)
        elif 'timeunit' in groupname2group:                 # per 1 time unit is implied as TimeInterval.
            quant = Quant(num_type = 'int', value = 1)
            time_unit_start = match_object.start('timeunit')
            time_unit = parse.position2struc(time_unit_start)
            time_interval = TimeInterval(quant = quant, time_unit = time_unit)
        elif debug:
            raise Exception("Error in parsing max dose in raw sig: ->%s<-\n Time Interval indicator is missing in found expression ->%s<-" %
                            (parse.sig.raw_sig, match_object.group()))
        else:
            return None

        if 'directive' in groupname2group:
            directive_start = match_object.start('directive')
            directive = parse.position2struc(directive_start)
        elif 'directive1' in groupname2group:
            directive_start = match_object.start('directive1')
            directive = parse.position2struc(directive_start)
        else:
            directive = None


        struc = MaxDose(dose = dose, time_interval = time_interval, directive = directive)

        return [struc]

    rule = Rule_ExtractStrucs(  name = 'max dose',
                                search_patterns = [pattern0, pattern1],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)
    return rule



def rule_assembly_stop_condition():

    pattern0 = re.compile(r'''
                            (?<![a-z])
                            (?P<until_gone>
                                (                       # pattern A: until gone/finished
                                    (
                                        untill?             |
                                        till?
                                    )
                                    \s
                                    (it \s is \s)?
                                    (all \s (are \s)? )?
                                    (
                                        gone                |
                                        finish(ed)?         |
                                        end                 |
                                        taken
                                    )
                                )                   |
                                (                       # pattern B: and finish all medicine
                                    ( (and|AND_CONJ) \s)?
                                    finish \s all
                                    (\s the)?
                                    (\s medicine)?
                                )
                            )
                            (?![a-z])
                        ''', re.X)

    pattern1 = re.compile(r'''
                            (?<![a-z])
                            (?P<until_completed>
                                (                               # pattern A: "till completion"
                                    (
                                        untill?     |
                                        till?
                                    )
                                    \s
                                    (
                                        complete    |
                                        completed   |
                                        completion
                                    )
                                )                           |
                                (                               # pattern B: "to completion" (we want to avoid "to complete"
                                    to \s completion
                                )
                            )
                            (?![a-z])
                        ''', re.X)

    pattern2 = re.compile(r'''                      # "To effect" (i.e. until the dose is effective).
                            (?<![a-z])
                            (?P<to_effect>
                                (
                                    (to|untill?|till?) \s
                                    (reach \s)?
                                    (maximal \s)?
                                    effect (ive|s)?             # mandatory in this option
                                    (\s dose)?
                                )                       |
                                (                               # from dictionary
                                    until \s the \s dose \s is \s effective
                                )
                            )
                            (?![a-z])
                        ''', re.X)






    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        groupname2group = trim_dictionary(match_object.groupdict())
        found = match_object.group()

        if 'until_gone' in groupname2group:
            condition = 'until_gone'
        elif 'until_completed' in groupname2group:
            condition = 'until_therapy_completed'
        elif 'to_effect' in groupname2group:
            condition = 'to_effect'
        else:
            return None

        struc = StopCondition(condition = condition, constituents = [found])

        return [struc]

    rule = Rule_ExtractStrucs(  name = 'stop_condition',
                                search_patterns = [pattern0, pattern1, pattern2],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)
    return rule


def rule_discard_remainder():
    pattern = re.compile(r'''
                            (?<![a-z])
                            (and \s|AND_CONJ \s)?
                            (then \s|THEN_CHRONO \s)?
                            discard \s                      # mandatory
                            (any \s|the \s)?
                            (
                                remainder                                                               |
                                (remaining|rmeaining|reamining|rmeaining|unused|remain)                 |
                                excess                                                                  |
                                extra                                                                   |
                                the \s rest
                            )
                            (\s portion| \s medication)?
                            (?![a-z])
                        ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        found = match_object.group()
        struc = DiscardRemanider(constituents = [found])

        return [struc]

    rule = Rule_ExtractStrucs(  name = 'discard_remainder',
                                search_patterns = [pattern],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)

    return rule


def rule_miscellanous_call_911():
    pattern0 = re.compile(r'''
                            (?<![a-z])
                            (?P<condition_before>
                                (   if \s not \s resolved:?                                         |
                                    if \s no \s relief                                              |
                                    if (\s chest)? \s pain \s (persists | cont(inues)?)             |
                                    if \s still \s chest \s pain (\s (persists | cont(inues)?))?    |
                                    if \s (?P<indication>INDICATION) \s (persists | cont(inues)?)   |
                                    (if \s)? doesn't \s help                                        |
                                    (if \s)? no \s help                                             |
                                    else
                                ) \,? \s?
                            )?
                            call \s?                # mandatory
                            (
                                (?P<num_911>QUANT)      |
                                (m\.?d\.? | doctor)
                            )
                            (?P<condition_after>
                                \s?
                                (   if \s not \s resolved:?                                         |
                                    if \s no \s relief                                              |
                                    if (\s chest)? \s pain \s (persists | cont(inues)?)             |
                                    if \s still \s chest \s pain (\s (persists | cont(inues)?))?    |
                                    if \s (?P<indication2>INDICATION) \s (persists | cont(inues)?)  |
                                    (if \s)? doesn't \s help                                        |
                                    (if \s)? no \s help
                                )
                            )?
                            (?![a-z])
                            ''', re.X)

    pattern1 = re.compile(r'''          # Remove redundant Then_Chrono which confuses the Sem parsing. E.g. "then if no relief call 911"
                            THEN_CHRONO \s (?P<misc>MISCELLANEOUS)
                            ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())

        if 'num_911' in groupname2group:
            num_911_start = match_object.start('num_911')
            num_911 = parse.position2struc(num_911_start)
            if num_911.value != 911:
                return
        if 'misc' in groupname2group:
            misc_start = match_object.start('misc')
            misc = parse.position2struc(misc_start)
            if misc.typ in ('call_911_if_condition_persists', 'call_911'):
                return [misc]
            else:
                return None

        if 'indication' in groupname2group:
            indication_start = match_object.start('indication')
            indication = parse.position2struc(indication_start)
            if indication.condition != 'pain':
                return
        elif 'indication2' in groupname2group:
            indication_start = match_object.start('indication2')
            indication = parse.position2struc(indication_start)
            if indication.condition != 'pain':
                return
        if 'condition_before' in groupname2group or 'condition_after' in groupname2group:
            typ = 'call_911_if_condition_persists'
        else:
            typ = 'call_911'

        found = match_object.group()
        struc = Miscellaneous(typ = typ, constituents = [found])

        return [struc]

    rule = Rule_ExtractStrucs(  name = 'Miscellaneous: Call 911',
                                search_patterns = [pattern0, pattern1],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)

    return rule

def rule_miscellanous_antibiotic():
    pattern0 = re.compile(r'''
                            (?<![a-z])          # should have either the prefix, e.g. "for antibiotic" or "(antibiotic)", etc.
                                                # or be at the end of the sig but NOT preceded by "after antibiotic" or "last day of antibiotic"e
                            (
                                (?P<prefix>
                                    is (\s an?)?        |
                                    as (\s an?)?        |
                                    for (\s an?)?       |
                                    \-                  |
                                    \*                  |
                                    \(
                                ) \s?
                            )?
                            anti (\s|\-)? bioti(c|k)s?   # mandatory
                            \*?
                            \)?
                            (?![a-z])
                            ''', re.X)


    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())
        if 'prefix' not in groupname2group:
            bad_left_context_pattern = re.compile(r'(after|of) \s$')
            if bad_left_context_pattern.search(left_context):
                return None
            elif trim(right_context) not in ('', ' ', '.'):
                return


        found = match_object.group()
        struc = Miscellaneous(typ = 'antibiotic', constituents = [found])

        return [struc]

    rule = Rule_ExtractStrucs(  name = 'Miscellaneous: antibiotic',
                                search_patterns = [pattern0],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)

    return rule


def rule_miscellanous_vitamin():
    pattern0 = re.compile(r'''
                            (?<![a-z])          # should not have a bad prefix, e.g. "take with vitamin a" or "take prior to vitamins", etc.
                                                # Always should be at the end of the sig. In almost all cases we saw, it is just one word "vitamin" out of any context at the very end.
                            (
                                (?P<bad_prefix>
                                    with            |
                                    prior \s to
                                ) \s?
                            )?
                            (
                                (?P<good_prefix>
                                    is (\s an?)?    |
                                    for
                                ) \s?
                            )?
                            (multi \-? \s?)?
                            vitamins?                   # mandatory
                            (\s [a-z])?                 # E.g. "vitamin a"
                            (?![a-z])
                            ''', re.X)


    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())
        if 'bad_prefix' in groupname2group:
            return None
        if trim(right_context) not in ('', ' ', '.'):
            return
        found = match_object.group()
        struc = Miscellaneous(typ = 'vitamin', constituents = [found])

        return [struc]

    rule = Rule_ExtractStrucs(  name = 'Miscellaneous: vitamin',
                                search_patterns = [pattern0],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)

    return rule

def rule_miscellanous_supplement():
    pattern0 = re.compile(r'''
                            (?<![a-z])          #
                                                # Always should be at the end of the sig. In almost all cases we saw, it is just one word "supplement" out of any context at the very end.
                            (
                                (?P<good_prefix>
                                    is (\s an?)?    |
                                    for
                                ) \s?
                            )?
                            (dietary  \s?)?
                            supplement s?                   # mandatory
                            (?![a-z])
                            ''', re.X)


    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())
        if trim(right_context) not in ('', ' ', '.'):
            return
        found = match_object.group()
        struc = Miscellaneous(typ = 'supplement', constituents = [found])

        return [struc]

    rule = Rule_ExtractStrucs(  name = 'Miscellaneous: supplement',
                                search_patterns = [pattern0],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)

    return rule



def rule_miscellanous_diuretic():
    pattern0 = re.compile(r'''
                            (?<![a-z])          # should not have a bad prefix, e.g. "take with water"
                                                # Always should be at the end of the sig. In almost all cases we saw, it is just one word "water" or "water pill" out of any context at the very end.
                            (
                                (?P<bad_prefix>
                                    with    |           # e.g. "take with water" or "take with a glass of water", "other than water"
                                    of      |
                                    than
                                )
                                (\s plain)?
                                \s?
                            )?
                            (
                                (?P<good_prefix>
                                    is (\s an?)?
                                ) \s?
                            )?
                            (                           # mandatory
                                diuretic s?         |
                                water (\s pills?)?
                            )
                            (?![a-z])
                            ''', re.X)


    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())
        if 'bad_prefix' in groupname2group:
            return None
        if trim(right_context) not in ('', ' ', '.'):
            return
        found = match_object.group()
        struc = Miscellaneous(typ = 'diuretic', constituents = [found])

        return [struc]

    rule = Rule_ExtractStrucs(  name = 'Miscellaneous: diuretic',
                                search_patterns = [pattern0],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)

    return rule




def rule_miscellanous_insulin():
    pattern0 = re.compile(r'''
                            (?<![a-z])          # "injections" is after "insulin": e.g "for lantus or humalog injections"
                            (?P<insulin_injection>
                                (
                                    (with|for|in|via)       # e.g. "in insulin pump" or "for insulin"
                                    \s
                                )
                                (
                                    humalog                 |
                                    lantus (\s solostar)?   |
                                    levemir                 |
                                    novolin                 |
                                    novolog                 |
                                    victoza                 |
                                    insulin s?              |
                                    insukin
                                )
                                (\s (and|AND_CONJ|or) \s (lantus|levemir|humalog|novolin|novolog|victoza) )?
                                (\s insulin)?
                                (?P<injections>             # "injection" is mandatory unless end of sig (e.g. "use as directed with lantus") or preceded by "use" or "use as directed"
                                    \s
                                    (
                                        injections?         |
                                        injectios?          |
                                        (flex)? pens?       |
                                        cartridges?         |
                                        pump
                                    )
                                )?
                            )
                            (?![a-z])
                            ''', re.X)

    pattern1 = re.compile(r'''
                            (?<![a-z])          # "inject" before "insulin". E.g. "to inject insulin"
                            (?P<inject_insulin>
                                (
                                    to
                                    \s
                                )?
                                (                                   # mandatory
                                    (?P<inject>DIRECTIVE)   |
                                    for \s injection \s     |
                                    inject
                                )
                                \s
                                (
                                    humalog                 |
                                    lantus (\s solostar)?   |
                                    levemir                 |
                                    novolin                 |
                                    novolog                 |
                                    victoza                 |
                                    insulin s?              |
                                    insukin
                                )
                                (\s insulin s?)?
                            )
                            (?![a-z])
                            ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())

        found = match_object.group()
        struc = Miscellaneous(typ = 'for_insulin_injections', constituents = [found])

        if 'insulin_injection' in groupname2group:                  # pattern0
            if 'injections' not in groupname2group and right_context not in ('', ' ', '.'):
                # If "injections" is not present in the phrase, make sure that the phrase is terminal (e.g. "use once daily with lantus."
                # or it is immediately preceded by "use" directive or AS_DIRECTED.
                good_left_context_pattern = re.compile(r'((?P<directive>DIRECTIVE)|(?P<as_directed>AS_DIRECTED)|use) \s$')
                left_obj = good_left_context_pattern.search(left_context)
                if not left_obj:
                    return None
                left_groupname2group = trim_dictionary(left_obj.groupdict())
                if 'directive' in left_groupname2group:
                    left_directive_start = left_obj.start('directive')
                    left_directive = parse.position2struc(left_directive_start)
                    if left_directive.value != 'use':
                        return None
                elif 'as_directed' in left_groupname2group:
                    pass
        elif 'inject_insulin' in groupname2group:                   # pattern1
            if 'inject' in groupname2group:
                directive_start = match_object.start('inject')
                directive = parse.position2struc(directive_start)
                if directive.value != 'inject':
                    return None
                else:                                       # return back directive, too. e.g. "inject insulin twice daily"
                    return [struc, Struc(' ', []), directive]




        return [struc]

    rule = Rule_ExtractStrucs(  name = 'Miscellaneous: insulin',
                                search_patterns = [pattern0, pattern1],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)

    return rule



def rule_miscellanous_do_not_swallow():
    """  do_not_swallow
    """
    pattern = re.compile(r'''
                            (?<![a-z])
                            (?P<do_not_swallow>
                                (do \s not| don\'t)
                                \s
                                (
                                    (?P<swallow>DIRECTIVE)  |
                                    swallow
                                )
                            )
                            (?![a-z])
                            ''', re.X)


    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())
        if 'do_not_swallow' in groupname2group:
            if 'swallow' in groupname2group:
                directive_start = match_object.start('swallow')
                directive = parse.position2struc(directive_start)
                if directive.value != 'swallow':
                    return None
            misc_type = 'do_not_swallow'

        found = match_object.group()
        struc = Miscellaneous(typ = misc_type, constituents = [found])

        return [struc]

    rule = Rule_ExtractStrucs(  name = 'Miscellaneous: do not swallow',
                                search_patterns = [pattern],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)

    return rule


def rule_miscellanous_rinse_mouth_after_use():
    """   rinse_mouth_after_use
    """

    pattern = re.compile(r'''
                            (?<![a-z])          #  *rinse mouth with water after use*
                            (?P<rinse_mouth_after_use>
                                (                               # mandatory
                                    (?P<rinse>DIRECTIVE)    |
                                    rinse                   |
                                    wash (\s your)?
                                )
                                (?P<mouth>                  # if "mouth" is not present, we will check if the context is right for "mouth rinse" as opposed to rinsing off e.g. skin
                                    (\s the)? \s mouth      |
                                    \s (?P<site_mouth>SITE)
                                )?
                                (?P<water>                              # optional.  "with water" could  be a Substrate
                                    \s with \s water                |
                                    \s? (?P<substrate>SUBSTRATE)
                                )?
                                (?P<after>\s after)?
                                (\s each)? (\s every)?
                                (?P<use>
                                    \s (?P<use_dir>DIRECTIVE)       |
                                    \s use                          |
                                    \s using
                                )?
                                \** \.?
                            )
                            (?![a-z])
                            ''', re.X)


    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())

        if 'rinse' in groupname2group:
            rinse_start = match_object.start('rinse')
            rinse = parse.position2struc(rinse_start)
            if rinse.value != 'rinse':
                return None
        if 'use_dir' in groupname2group:
            use_start = match_object.start('use_dir')
            use = parse.position2struc(use_start)
            if use.value != 'use':
                return None
        if 'substrate' in groupname2group:
            substrate_start = match_object.start('substrate')
            substrate = parse.position2struc(substrate_start)
            if substrate.value != 'water':
                return None

        if 'site_mouth' in groupname2group:
            site_start = match_object.start('site_mouth')
            site = parse.position2struc(site_start)
            if site.value != 'mouth':
                return None

        left_pattern = re.compile('(?P<directive>DIRECTIVE)')
        left_obj = left_pattern.search(left_context)
        if 'mouth' not in groupname2group:
            # if "mouth" is not present, we will check if the context is right for "mouth rinse" as opposed to rinsing off e.g. skin
            # The left context is right if it has 'inhale'
            if not left_obj:
                return None
            left_groupname2group = trim_dictionary(left_obj.groupdict())
            if 'directive' in left_groupname2group:
                left_directive_start = left_obj.start('directive')
                left_directive = parse.position2struc(left_directive_start)
                if left_directive.value != 'inhale':
                    return None
        if 'water' not in groupname2group and 'use' not in groupname2group and 'after' not in groupname2group and 'after' not in right_context and 'after' not in left_context and 'THEN_CHRONO' not in left_context:
            # Eliminate the possibility of "rinse mouth twice a day for 30 seconds", which does NOT mean rinse after use.
            # Or should have "after" somewhere or be "then rinse". E.g., avoid "use as directed to rinse mouth"
            # In this case make sure that there is some other directive on the left, or quit
            if not left_obj:
                return None

        misc_type = 'rinse_mouth_after_use'
        found = match_object.group()
        struc = Miscellaneous(typ = misc_type, constituents = [found])

        return [struc]

    rule = Rule_ExtractStrucs(  name = 'Miscellaneous: rinse mouth after use.',
                                search_patterns = [pattern],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)

    return rule

def rule_miscellanous_administer_in_office_or_pharmacy():
    pattern0 = re.compile(r'''                  # to be administed in doctor's office or pharmacy
                            (?<![a-z])
                                ((to \s be | will \s be | for ) \s)?
                                (administered | administration | given | used)
                                \s
                                (by|in|at)
                                \s
                                (the \s)? (your \s)?
                                (
                                    (physician | m\.?d\.? | doctor) (\'?s)? (\s office)?    |
                                    office                                                  |
                                    nurses? (\s practitioner)? (\s practioner)?             |
                                    (pharmacy | pharmacists?)
                                )
                            (?![a-z])
                            ''', re.X)

    pattern1 = re.compile(r'''                  # for office use
                            (?<![a-z])
                                (
                                    (?P<directive>DIRECTIVE)    # Take
                                    \s to \s
                                    (the \s)?
                                    (
                                        clinic                                          |
                                        (doctor | physician |md) \'?s? \s office        |
                                        (pharmacy | pharmacists?)
                                    )
                                \s)?
                                (for \s)?
                                (office | pharmacy | pharmacists?)     # mandatory
                                \s
                                (           # mandatory: USE  or administration
                                    (?P<use>DIRECTIVE)      |
                                    administration
                                )
                            (?![a-z])
                            ''', re.X)

    pattern2 = re.compile(r'''                  # use at/take to doctor's office
                            (?<![a-z])
                                (                                   # Optional: Take or Use or bring
                                    (?P<directive>DIRECTIVE)    |
                                    bring
                                \s)?
                                (at|to) \s                          # Mandatory: at (doctor|md etc) office
                                (the \s)?(your \s)?
                                (
                                    clinic                                              |
                                    (doctor | physician |md) \'?s? \s office            |
                                    office                                              |
                                    (pharmacy | pharmacists?)
                                )
                            (?![a-z])
                            ''', re.X)

    pattern3 = re.compile(r'''                  # misc: for in office use/sedation,
                            (?<![a-z])
                                (for \s)?
                                (the \s)?
                                (
                                    in (\s|\s?\-\s?) office     |
                                    pharmacy                    |
                                    pharmacists?
                                )
                                (\s
                                    (sedation | administration | admin | (?P<use>DIRECTIVE))
                                )
                            (?![a-z])
                            ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())

        if 'directive' in groupname2group:
            directive_start = match_object.start('directive')
            directive = parse.position2struc(directive_start)
            if directive.value not in ('take', 'use'):
                return None
        if 'use' in groupname2group:
            use_start = match_object.start('use')
            use = parse.position2struc(use_start)
            if use.value != 'use':
                return None

        found = match_object.group()
        if 'pharm' in found:
            typ = 'administer_in_pharmacy'
        else:
            typ = 'administer_in_office'
        struc = Miscellaneous(typ = typ, constituents = [found])

        return [struc]

    rule = Rule_ExtractStrucs(  name = 'Miscellaneous: Administer in office',
                                search_patterns = [pattern0, pattern1, pattern2, pattern3],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)

    return rule

def rule_miscellanous_inject_in_office_or_pharmacy():
    pattern0 = re.compile(r'''                  # to be injected by doctor or pharmacist
                            (?<![a-z])
                                (to \s be | will \s be | for)
                                \s
                                (injected | injection |inject)
                                \s
                                (intramuscullarly \s)?
                                (by|in|at)
                                \s
                                (the \s)?
                                (
                                    (physician | m\.?d\.? | doctor) (\'s)? (\s office)?     |
                                    (pharmacy | pharmacists?)                               |
                                    office
                                )
                            (?![a-z])
                            ''', re.X)

    pattern1 = re.compile(r'''                  # if part of the pattern was already pre-processed by MISC of typ = administer_in_office
                            (?<![a-z])
                                (ANAPHORA \s)?
                                (to \s be | will \s be | for (\s the)? )
                                \s
                                (injected | injection)
                                \s
                                (intramuscullarly \s)?
                                (?P<prior_misc>MISCELLANEOUS)
                            (?![a-z])
                            ''', re.X)


    pattern2 = re.compile(r'''                  # for office injection
                            (?<![a-z])
                                (
                                    ((?P<directive>DIRECTIVE)|bring)    # Take
                                    \s to \s
                                    (the \s)?
                                    (
                                        clinic                                          |
                                        (doctor | physician |md) \'?s? \s office        |
                                        (pharmacy | pharmacists?)
                                    )
                                \s)?
                                for \s              # mandatory
                                (the \s)?
                                (office | pharmacist | pharmacy)    # mandatory
                                \s?
                                injection           # mandatory
                            (?![a-z])
                            ''', re.X)

    pattern3 = re.compile(r'''                  # to inject in doctor's office
                            (?<![a-z])
                                ( (to|for) \s)?
                                (                   # Mandatory: Inject
                                    (
                                        (?P<inject>DIRECTIVE)   |
                                        inject                  |
                                        shot                                # yes, there is an expression "shot in pharmacy"
                                    )
                                    \s
                                )
                                (at|in|by) \s
                                (the \s)?
                                (
                                    clinic                                          |
                                    ((doctor | physician |md) \'?s? \s)? office     |
                                    (pharmacy | pharmacists?)
                                )
                            (?![a-z])
                            ''', re.X)
    pattern4 = re.compile(r'''                  # "in office" if "inject" is in the left context or in the INJECT AS DIRECED
                            (?<![a-z])
                                (?P<inject_in_left_context>
                                    in \s (the \s)? office              |
                                    (
                                        ((at|in|by|for) \s )
                                        (the \s)?
                                        (pharmacy | pharmacists?)
                                    )
                                )
                            (?![a-z])
                            ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        groupname2group = trim_dictionary(match_object.groupdict())

        prior_misc = None

        if 'inject' in groupname2group:
            inject_start = match_object.start('inject')
            inject = parse.position2struc(inject_start)
            if inject.value != 'inject':
                return None

        if 'directive' in groupname2group:
            directive_start = match_object.start('directive')
            directive = parse.position2struc(directive_start)
            if directive.value not in ('take', 'use'):
                return None
        if 'prior_misc' in groupname2group:       # pattern1
            misc_start = match_object.start('prior_misc')
            prior_misc = parse.position2struc(misc_start)
            if prior_misc.typ not in ('administer_in_office', 'inject_in_office', 'administer_in_pharmacy', 'inject_at_pharmacy'):
                return None

        if 'inject_in_left_context' in groupname2group:     # pattern4
            left_pattern = re.compile('(?P<directive>DIRECTIVE)|(?P<as_directed>AS_DIRECTED)')
            left_obj = left_pattern.search(left_context)
            if not left_obj:
                return None
            left_groupname2group = trim_dictionary(left_obj.groupdict())
            if 'directive' in left_groupname2group:
                left_directive_start = left_obj.start('directive')
                left_directive = parse.position2struc(left_directive_start)
                if left_directive.value != 'inject':
                    return None
            else:                   # AS_DIRECTED is in the left context
                left_as_directed_start = left_obj.start('as_directed')
                left_as_directed = parse.position2struc(left_as_directed_start)
                if left_as_directed.directive != 'inject':
                    return None

        found = match_object.group()
        if 'pharm' in found or prior_misc and 'pharm' in prior_misc.typ:
            typ = 'inject_at_pharmacy'
        else:
            typ = 'inject_in_office'
        struc = Miscellaneous(typ = typ, constituents = [found])

        return [struc]

    rule = Rule_ExtractStrucs(  name = 'Miscellaneous: Inject in office',
                                search_patterns = [pattern0, pattern1, pattern2, pattern3, pattern4],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)

    return rule



def rule_and_conj():
    pattern = re.compile(r'''
                            (?<![a-z])
                            (
                                (and \s)? also    |
                                and
                            )
                            (?P<spit>       # ... "and spit out" is part of the previous event and should either be merged with previous directive ("swish and spit") or just deleted out.
                                \s
                                (
                                    spit (\s out)?   |
                                    expectorate
                                )
                            )?
                            (?![a-z])
                        ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        groupname2group = trim_dictionary(match_object.groupdict())

        found = match_object.group()
        struc = AndConj(constituents = [found])

        if 'spit' in groupname2group:
            # ... "and spit out" is part of the previous event and should either be merged with previous directive ("swish and spit") or just deleted out.
            # E.g. "rinse half-ounce in mouth for 30 seconds and spit out twice a day"
            left_context_pattern = re.compile('DIRECTIVE')
            prior_directive_obj = left_context_pattern.search(left_context)
            if not prior_directive_obj:
                return []
            prior_directive_start = prior_directive_obj.start()
            prior_directive = parse.position2struc(prior_directive_start)

            if prior_directive.value in ('swish', 'rinse'):
                prior_directive.value = 'swish and spit'
                prior_directive.rules_used.append('deduced_in_then_chrono_rule')
            return []

        return [struc]

    rule = Rule_ExtractStrucs(  name = 'and_conjunction',
                                search_patterns = [pattern],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)

    return rule


def rule_then_chrono():
    pattern = re.compile(r'''
                            (?<![a-z])
                            ((and|AND_CONJ) \s)?
                            (?P<to_start>                   # "to start" can have 2 meanings: "to start 2 days before surgery" is a specification of timing. But "take 2 tabs to start then increase..." means "on day 1: take 2 tabs"
                                to \s start \s (with \s)?
                            )?
                            (
                                then,?          |
                                after \s that:?
                            )
                            (?P<spit>       # ... "then spit out" is part of the previous event and should either be merged with previous directive ("swish and spit") or just deleted out.
                                \s
                                (
                                    spit (\s out)?   |
                                    expectorate
                                )
                            )?
                            (?![a-z])
                        ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        groupname2group = trim_dictionary(match_object.groupdict())

        found = match_object.group()
        struc = ThenChrono(constituents = [found])

        if 'spit' in groupname2group:
            # ... "then spit out" is part of the previous event and should either be merged with previous directive ("swish and spit") or just deleted out.
            # E.g. "rinse half-ounce in mouth for 30 seconds then spit out twice a day"

            left_context_pattern = re.compile('DIRECTIVE')
            prior_directive_obj = left_context_pattern.search(left_context)
            if not prior_directive_obj:
                return []
            prior_directive_start = prior_directive_obj.start()
            prior_directive = parse.position2struc(prior_directive_start)

            if prior_directive.value in ('swish', 'rinse'):
                prior_directive.value = 'swish and spit'
                prior_directive.rules_used.append('deduced_in_then_chrono_rule')
            return []
        if 'to_start' in groupname2group:
            # "to start" can have 2 meanings: "to start 2 days before surgery" is a specification of timing. But "take 2 tabs to start then increase..." means "on day 1: take 2 tabs"
            cal_event = Calendar_Event(typ = 'relative')
            cal_event.value = 'now'
            cal_event.rules_used.append('*deduced_in_then_chrono_rule*')
            return [cal_event, Struc(' ', []), struc]

        return [struc]

    rule = Rule_ExtractStrucs(  name = 'then_chrono',
                                search_patterns = [pattern],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)

    return rule


def rule_anaphora():
    pattern = re.compile(r'''
                            (?<![a-z])
                            (
                                (?P<this_medicine>(of \s)? this \s (medication|medicine))   |
                                (?P<this_dose>this \s dose)                                 |
                                (?P<do_this>(please \s)? do \s this)                        |
                                (?P<this_is>this \s is)                                     |
                                (?P<this_will>this \s will)
                            )
                            (?![a-z])
                        ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        found = match_object.group()
        groupname2group = trim_dictionary(match_object.groupdict())
        value = groupname2group.keys()[0]
        struc = Anaphora(value = value, constituents = [found])

        return [struc]

    rule = Rule_ExtractStrucs(  name = 'anaphora',
                                search_patterns = [pattern],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)

    return rule



def rule_on_before_off():
    """ Cases such as "apply 12 hours on and 12 hours off daily 1 patch".  Here we deal just with the "on" part (assuming "off" is somewhere to the right.

    "12 hours on" really means DURATION, and "12 hours off" really means "Remove or stop for DURATION"

    So we need to insert correct DURATION.

    Also see process_special_duration_cases()
    """

    pattern0 = re.compile(r'''                      # e.g. "every 12 hours on" or "(for) 12 hours on"
                                (?P<prefix>DURATION|TIMEINTERVAL|PERIODICITY)
                                \s
                                (leave \s)?
                                on
                                (?![a-z])
                        ''', re.X)

    pattern1 = re.compile(r'''                      # or e.g. "on 12 hours"
                                (?<![a-z])
                                (                   # mandatory
                                    (leave \s)? on \s           |
                                    leave \s in \s place \s
                                )
                                (?P<suffix>DURATION|TIMEINTERVAL)
                        ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        correct_right_context = re.compile(r'''   # the right context has to have "off" or DIRECTIVE = "remove"  e.g. "then off 12 hours" or "and 12 hours off" or "then remove for 12 hours"
                                    (
                                        ((?P<right_duration_0>DURATION)|(?P<right_time_interval_0>TIMEINTERVAL)|(?P<right_periodicity>PERIODICITY)) \s (?P<directive_remove_0>DIRECTIVE)
                                    |
                                        (?P<directive_remove_1>DIRECTIVE)       # could be just "off" as in "12 weeks on, off for 12 weeks"
                                        \s
                                        ((?P<right_duration_1>DURATION)|(?P<right_time_interval_1>TIMEINTERVAL))
                                    )
                            ''', re.X)

        found_right_context = correct_right_context.search(right_context)
        if not found_right_context:
            return None
        else:
            right_groupname2group = trim_dictionary(found_right_context.groupdict())
            for i in range(2):      # check that if groups directive_remove_0, directive_remove_1 capture anything that thing is directive "remove"
                group_name = 'directive_remove_'  + str(i)
                if group_name in right_groupname2group:
                    directive_start = found_right_context.start(group_name) + len(left_context) + len(match_object.group())
                    directive_remove = parse.position2struc(directive_start)
                    if not directive_remove.value in ('remove', 'stop'):
                        return None

        groupname2group = trim_dictionary(match_object.groupdict())

        if 'prefix' in groupname2group:
            struc_start = match_object.start('prefix')
        else:
            struc_start = match_object.start('suffix')
        struc = parse.position2struc(struc_start)
        if struc.label == 'DURATION':
            duration = struc
        elif struc.label == 'TIMEINTERVAL':
            time_interval = struc
            duration = Duration(time_interval)
        elif struc.label == 'PERIODICITY':
            # Check that we are not in something like "apply every day on for 12 hours...", i.e. that the Periodicity before "on" refers to the plausible time, not to "daily"
            if 'prefix' in groupname2group:
                if 'right_duration_0' in right_groupname2group:
                    right_duration_start = found_right_context.start('right_duration_0') + len(left_context) + len(match_object.group())
                    right_duration = parse.position2struc(right_duration_start)
                    right_time_unit_value = right_duration.time_interval.time_unit.value
                elif 'right_duration_1' in right_groupname2group:
                    right_duration_start = found_right_context.start('right_duration_1') + len(left_context) + len(match_object.group())
                    right_duration = parse.position2struc(right_duration_start)
                    right_time_unit_value = right_duration.time_interval.time_unit.value
                elif 'right_time_interval_0' in right_groupname2group:
                    right_time_interval_start = found_right_context.start('right_time_interval_0') + len(left_context) + len(match_object.group())
                    right_time_interval = parse.position2struc(right_time_interval_start)
                    right_time_unit_value = right_time_interval.time_unit.value
                elif 'right_time_interval_1' in right_groupname2group:
                    right_time_interval_start = found_right_context.start('right_time_interval_1') + len(left_context) + len(match_object.group())
                    right_time_interval = parse.position2struc(right_time_interval_start)
                    right_time_unit_value = right_time_interval.time_unit.value
                elif 'right_periodicity' in right_groupname2group:
                    right_periodicity_start = found_right_context.start('right_periodicity') + len(left_context) + len(match_object.group())
                    right_periodicity = parse.position2struc(right_periodicity_start)
                    right_time_unit_value = right_periodicity.time_interval.time_unit.value
                else:
                    right_time_unit_value = None
                if right_time_unit_value and right_time_unit_value != struc.time_unit.value:
                    return None
            time_interval = TimeInterval(struc.quant, struc.time_unit)
            duration = Duration(time_interval)
        else:
            return None

        if duration.time_interval.time_unit.value == 'hour':
            duration.on_off = 'on'

        return [duration]

    rule = Rule_ExtractStrucs(  name = 'on_before_off',
                                search_patterns = [pattern0, pattern1],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)

    return rule


def rule_off_after_on():
    """ Cases such as "apply 12 hours on and 12 hours off daily 1 patch".  Here we deal just with the "off" part (assuming "on" is somewhere to the left.

    "12 hours off" really means "Remove or stop DURATION"

    So we need to insert correct DURATION and if needed insert a THEN_CHRONO
    Also see process_special_duration_cases()
    """

    pattern0 = re.compile(r'''   # e.g. "and 12 hours off"
                                (
                                    ((?P<directive_remove_0>DIRECTIVE) \s FORM \s)?    # e.g. "remove patch then leave off for 12 hours"
                                    (
                                        AND_CONJ (\s THEN_CHRONO)?  |
                                        THEN_CHRONO
                                    ) \s
                                )?
                                (?P<prefix>DURATION|TIMEINTERVAL|PERIODICITY) \s (?P<directive_remove_1>DIRECTIVE)

                        ''', re.X)

    pattern1 = re.compile(r'''   # e.g. "then off 12 hours" or  "then remove for 12 hours"
                                (
                                    ((?P<directive_remove_0>DIRECTIVE) \s FORM \s)?    # e.g. "remove patch then leave off 12 hours"
                                    (
                                        AND_CONJ (\s THEN_CHRONO)?  |
                                        THEN_CHRONO
                                    ) \s
                                )?
                                (?P<directive_remove_2>DIRECTIVE)               # could be just "off" as in "12 weeks on, off for 12 weeks"
                                \s
                                (?P<suffix>DURATION|TIMEINTERVAL)
                        ''', re.X)

    pattern2 = re.compile(r'''   # e.g. "remove and allow 1 week rest"
                                (?P<directive_remove_2>DIRECTIVE)
                                ,? \s? (AND_CONJ \s)?
                                allow \s
                                (?P<suffix>DURATION|TIMEINTERVAL)
                                (\s rest)?
                        ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):

        # verify that the "off" part is preceded by a an "on" part. If DURATION is found on the left, but it is not "on",
        # all is not lost, because we have cases such as "APPLY 1 PATCH TOPICALLY FOR 12 HOURS.REMOVE PATCH AND LEAVE OFF FOR 12 HOURS"
        # where the "on" part doesn't say explicitly "on". So after we identified the timeinterval in the current Duration,
        # we go back to the preceding duration and, if they both are "12 hours", we then change the preceding duration to on_off = 'on"

        correct_left_context = 'DURATION'
        if not correct_left_context in left_context:
            return None
        groupname2group = trim_dictionary(match_object.groupdict())

        for i in range(3):      # check that if groups directive_remove_0, directive_remove_1, directive_remove_2 capture anything that thing is directive "remove"
            group_name = 'directive_remove_'  + str(i)
            if group_name in groupname2group:
                directive_start = match_object.start(group_name)
                directive_remove = parse.position2struc(directive_start)
                if not directive_remove.value in ('remove', 'stop'):
                    return None
        directive = directive_remove        # either directive_remove_1 or directive_remove_2 have to be present, so directive is one of these 2.

        if 'prefix' in groupname2group:
            struc_start = match_object.start('prefix')
        else:
            struc_start = match_object.start('suffix')
        struc = parse.position2struc(struc_start)
        if struc.label == 'DURATION':
            duration = struc
        elif struc.label == 'TIMEINTERVAL':
            time_interval = struc
            duration = Duration(time_interval)
        elif struc.label == 'PERIODICITY':
            time_interval = TimeInterval(struc.quant, struc.time_unit)
            duration = Duration(time_interval)
        else:
            return None

        # Find the right-most Duration struc to the left (there are cases of multiple Durations, e.g. "apply 1 patch for 30 days 12 hours on 12 hours off")
        duration_on_start = left_context.rfind(correct_left_context)
        duration_on = parse.position2struc(duration_on_start)
        if duration_on.on_off != 'on':
            # last possibility: they are both "12 hours", then we make the duration_on to be on_off = "on"
            duration_on_timeinterval = duration_on.time_interval
            current_timeinterval = duration.time_interval
            if duration_on_timeinterval.quant.value == 12 and duration_on_timeinterval.time_unit.value == 'hour' and  current_timeinterval.quant.value == 12 and current_timeinterval.time_unit.value == 'hour':
                duration_on.on_off = 'on'
            elif duration_on_timeinterval.quant.num_type == 'var' and duration_on_timeinterval.time_unit.value == 'hour' and  current_timeinterval.quant.num_type == 'var' and current_timeinterval.time_unit.value == 'hour':
                # from dictionary
                duration_on.on_off = 'on'

        if duration.time_interval.time_unit.value == 'hour':
            duration.on_off = 'off'

        # To avoid infinite loops, if we are already dealing with "Remove Duration", after adjusting Duration's on_off attribute and the preceding Duration's
        # on_off attribute, return
        if struc.label == 'DURATION':
            return

        then_struc = ThenChrono(constituents = [])
        space0 = Struc(label = ' ')
        space1 = Struc(label = ' ')

        return [then_struc, space0, directive, space1, duration]

    rule = Rule_ExtractStrucs(  name = 'off_after_on',
                                search_patterns = [pattern0, pattern1, pattern2],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)

    return rule

def rule_remove_parenthesis_duplicate_quant():
    """ Deal with cases such as "Take two (2) tabs..."
    """

    pattern = re.compile(r'''
                            (?P<quant1> QUANT)
                            \s?
                            \(
                            (?P<quant2> QUANT)
                            \s?
                            \)
                        ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        quant1_start = match_object.start('quant1')
        quant1 = parse.position2struc(quant1_start)

        quant2_start = match_object.start('quant2')
        quant2 = parse.position2struc(quant2_start)

        if quant1.value and quant1.value == quant2.value:
            return [quant1]
        else:
            return

    rule = Rule_ExtractStrucs(  name = 'remove_parenthesis_duplicate_quant',
                                search_patterns = [pattern],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)

    return rule


def rule_cleanup_ordinals_quant():
    """ Deal with cases such as "Take  on 1st thru 12th days of the month" where the ordinals are given in mixed form.

    This really should be dealt with in the numeric_id.py module, but the way it is currently used it is hard to do.
    """

    pattern = re.compile(r'''
                            (?P<quant>QUANT)
                            \s?
                            (?P<suffix>st|nd|rd|th)
                            (?![a-z])
                        ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        quant_start = match_object.start('quant')
        quant = parse.position2struc(quant_start)

        suffix = match_object.group('suffix')

        if quant.num_type != 'int':
            return None

        if quant.value  % 10 == 1 and quant.value % 100 != 11:  # "st" suffix only works for 1, 21, 31, .. 101, ... but not for numbers = 11 mod 100 (e.g. not for 11)
            if suffix == 'st':
                pass
            else:
                return None
        elif quant.value % 10 == 2 and quant.value % 100 != 12:
            if suffix == 'nd':
                pass
            else:
                return None
        elif quant.value % 10 == 3 and quant.value % 100 != 13:
            if suffix == 'rd':
                pass
            else:
                return None
        elif suffix == 'th':
            pass
        else:
            return None

        quant.num_type = 'ord'
        return [quant]


    rule = Rule_ExtractStrucs(  name = 'cleanup_ordinal_quants',
                                search_patterns = [pattern],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)

    return rule

def rule_clean_up_punctuation_1():
    pattern1 = re.compile(r'''                  # e.g. " (MAXDOSE)"
                            \(
                            \s?
                            (?P<struc> MAXDOSE | AS_DIRECTED | DISCARD_REMAINDER )
                            \s?
                            \)
                            ''', re.X)

    pattern2 = re.compile(r'''                  # Terminal " -- MAXDOSE"
                            (\-)+
                            \s?
                            (?P<struc> MAXDOSE | AS_DIRECTED | DISCARD_REMAINDER )
                            \s?
                            $
                            ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        struc_start = match_object.start('struc')
        struc = parse.position2struc(struc_start)
        if left_context and left_context[-1].isupper():
            space_struc = Struc(' ')
            return [space_struc, struc]
        else:
            return [struc]

    rule = Rule_ExtractStrucs(  name = 'clean_up_punctuation_1',
                                search_patterns = [pattern1, pattern2],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)

    return rule

def rule_clean_up_punctuation_2():
    """ remove terminal periods and multiple dashes     """

    pattern = re.compile(r'''
                            \. \s? $
                            |
                            (?P<dashes>\-\-+)
                        ''', re.X)

    def search_proc(txt, search_pattern, start = None):
        start_pos = start if start else 0
        match_obj = search_pattern.search(txt, start_pos)
        return match_obj

    def replacement_proc(match_object, left_context, right_context, parse):
        if 'dashes' in match_object.groupdict() and right_context and right_context[0] != ' ':
            space_struc = Struc(' ')
            return [space_struc]
        else:
            return []

    rule = Rule_ExtractStrucs(  name = 'clean_up_punctuation_2',
                                search_patterns = [pattern],
                                search_proc = search_proc,
                                replacement_proc = replacement_proc,
                                rule_properties = set(['struc_identification']),
                                change_parse_in_place = True)

    return rule



RULES = (   rule_numerics_identification(),
            rule_cleanup_ordinals_quant(),
            rule_remove_parenthesis_duplicate_quant(),
            rule_assembly_varquant(),
            rule_assembly_as_directed_1(),
            rule_assembly_as_directed_2(),
            rule_assembly_form(),
            rule_assembly_directive(),
            rule_assembly_broken_up_directives(),
            rule_assembly_dir_attribute(),
            rule_assembly_range_identifcation(),
            rule_assembly_dose(),
            rule_remove_parenthesis_duplicate_dose(),
            rule_assembly_substrate(),
            rule_assembly_time_unit(),
            rule_assembly_timeinterval(),
            rule_assembly_freq(),
            rule_assembly_periodicity(),
            rule_assembly_route(),
            rule_assembly_duration(),
            rule_assembly_repeat(),
            rule_assembly_timing(),
            rule_assembly_time_of_day(),
            rule_assembly_calendar_event(),            # rule_assembly_nightly(),
            rule_assembly_specific_day(),
            rule_assembly_indication_pain(),
            rule_assembly_indication(),
            rule_assembly_apply_to_site(),
            rule_assembly_instillables_site(),
            rule_assembly_other_site(),
            rule_assembly_as_needed(),
            rule_assembly_see_instructions(),
            rule_assembly_vehicle(),
            rule_assembly_dose_omitted_form(),
            rule_assembly_spray_n_times(),
            rule_assembly_maxdose(),
            rule_assembly_stop_condition(),
            rule_discard_remainder(),
            rule_assembly_taper(),
            rule_miscellanous_call_911(),
            rule_miscellanous_antibiotic(),
            rule_miscellanous_vitamin(),
            rule_miscellanous_supplement(),
            rule_miscellanous_diuretic(),
            rule_miscellanous_insulin(),
            rule_miscellanous_do_not_swallow(),
            rule_miscellanous_rinse_mouth_after_use(),
            rule_miscellanous_administer_in_office_or_pharmacy(),
            rule_miscellanous_inject_in_office_or_pharmacy(),
            rule_and_conj(),
            rule_then_chrono(),
            rule_anaphora(),
            rule_on_before_off(),
            rule_off_after_on(),
            rule_clean_up_punctuation_1(),
            rule_clean_up_punctuation_2()
        )



def cleanup_raw_sig(raw_sig):
    """ Clean raw sig by normalizing and doing spelling corrections. Then create and return sig object. """

    raw_sig = str(raw_sig)
    cleaned_raw_sig = normalize_string(raw_sig)
    cleaned_raw_sig = spelling_corrections(cleaned_raw_sig)
    cleaned_raw_sig = trim(cleaned_raw_sig)
    cleaned_raw_sig = latin_2_sig(cleaned_raw_sig)
    cleaned_raw_sig = trim(cleaned_raw_sig)

    sig = Sig(raw_sig, cleaned_raw_sig)

    return sig

def apply_preprocess_rules(sig):
    """ PREPROCESS
        A small group of "local",  context-independent Rules that need to apply before the larger patterns can be determined.
        Includes extraction of numerical expressions, abbreviations, spelling normalization.
    """

    preprocessing_rules = [rule for rule in RULES if 'preprocessing' in rule.rule_properties]
    for parse in sig.parses:
        for rule in preprocessing_rules:
            rule.apply_rule(sig)


def apply_struc_identification_rules(sig):
    """ Apply main rules to create Strucs """

    struc_identification_rules = [rule for rule in RULES if 'struc_identification' in rule.rule_properties]

    def one_run_through_all_rules():
        for rule in struc_identification_rules:
            rule.apply_rule(sig)

    for parse in sig.parses:
        parse.changed_on_last_pass = False
        for run_num in range(10):       # just for safety reasons we limit to 10 passes, in case the changes are circular.
            parse.changed_on_last_pass = False
            one_run_through_all_rules()
            parse.number_of_cycles_of_rule_application += 1
            if not parse.changed_on_last_pass:
                break



    if debug:
        for parse in sig.parses:
            for struc in parse.strucs:
                valid = struc.is_valid_struc()
                if not valid:
                    raise Exception('for sig -->%s<-- struc not valid: \n%s' % (sig.raw_sig, struc))
        pass
        #parse = sig.parses[0]
        #strucs = parse.strucs
        #print('\n------------\n Raw Sig: %s' % sig.raw_sig)
        #start = 0
        #for i, struc in enumerate(strucs):
        #    end = start + len(struc.label)
        #    print('Struc %2d. Start: %3d  End: %3d Label: %-20s  Value: |%s|' % (i, start, end, struc.label, struc.value))
        #    start = end
        #print('\n------------\n Flattened Parse: Len: %d, Parse: \n%s' % (len(parse.flatten()), parse.flatten()))



class SemScopeObject(object):
    """ An informational structure assigned to some Strucs (those that are crucial to deciding if we need to start a new Schedule or a new Event).
    """
    def __init__(self, num, label):
        self.struc_index = num
        self.struc_label = label
        self.is_duration_like = None
        self.is_timing_like = None
        self.left_scope_start = None
        self.right_scope_end = None
        self.is_right_scoped = None
        self.is_left_scoped = None

    def pprint(self):
        pairs = sorted([prop + ':' + str(val) for prop, val in self.__dict__.items() if val is not None], reverse = True)
        return '|'.join(pairs)


class SemScopeData(object):
    """ Provides information on the scope of semantic objects such as Instruction, Schedule, and AdminEvent.

        We need to know for each struc if it should it start a new Schedule, AdminEvent, or nothing.

        The struc needs to start a new Schedule if:
        - If the struc IS right-scoped Duration or Then_Chrono immediately preceding the right_scoped Duration.
          Also then assign a scope of that Duration until the next Then_Chrono or Duration.
        - Elif the struc is not in the scope of a right-scoped Duration and there is a Duration sctruc on the right which has no scope start yet.
          Also then assign the struc to be the start of the scope of that Duration struc.
        - Elif There is no current_schedule.

        The struc needs to start a new AdminEvent if:
        - We are forced to start something new (Schedule or AdminEvent) because struc == Then_Chrono or an incompatible struc
            AND we are within the scope of an existing Schedule and don't need to start a new one
            AND there is a Timing or SPECIFIC_DAY on the right that doesn't have a scope start yet.
          In this case, start AdminEvent and assign the scope to the Timing
        - Elif we are an AdminEvent-type struc and there is no current_event (== None)

    Note on 4/26/2013: The property is_right_scoped seems to be unnecessary and is never used. All duration-like strucs seem to be left scoped. But it works in gold_standard cases and
    I don't have the time to fix the code now.
    Where it doesn't work and we need to fix the right_scoped rule is a case such as
    "Take 2 tabs in AM and 3 in PM for 4 days. Then for the next 7 days: take 1 tab in the AM and 1.5 in the PM. Then for the next 30 days: take 6 tab daily."
    LF

    """

    def __init__(self, struc_list):
        self.duration_like_strucs = []
        self.timing_like_strucs = []

        for (struc_index, struc) in enumerate(struc_list):
            if struc.label in ('DURATION', 'CALENDAR_EVENT', 'TAPER'):
                self.close_scope_for_right_scoped_duration(struc_index)
                sem_scope_object = SemScopeObject(struc_index, struc.label)
                sem_scope_object.is_duration_like = True
                self.duration_like_strucs.append(sem_scope_object)
                if struc_index > 1 and struc_list[struc_index - 2] == 'THEN_CHRONO':
                    # Note on 4/26/2013: The property is_right_scoped seems to be unnecessary and is never used. The above condition is never true because .label is missing.
                    # But the algorithm works fine without it. And I don't have time to fix it now. LF
                    sem_scope_object.is_right_scoped = True
                elif struc_index + 1 < len(struc_list) and struc_list[struc_index + 1].label == ':':
                    sem_scope_object.is_right_scoped = True
                else:
                    sem_scope_object.is_left_scoped = True

            elif struc.label == 'THEN_CHRONO':
                self.close_scope_for_right_scoped_duration(struc_index)
            elif struc.label in ('TIMING', 'SPECIFIC_DAY'):
                sem_scope_object = SemScopeObject(struc_index, struc.label)
                sem_scope_object.is_timing_like = True
                self.timing_like_strucs.append(sem_scope_object)

    def close_scope_for_right_scoped_duration(self, num):
        if self.duration_like_strucs and self.duration_like_strucs[-1].is_right_scoped:
            self.duration_like_strucs[-1].right_scope_end = num - 1

    def get_list_member_with_largest_index_less_than_num(self, num, alist):
        last_object = None
        for sem_scope_object in alist:
            if num >= sem_scope_object.struc_index:
                last_object = sem_scope_object
            else:
                return last_object
        # in case alist is empty:
        return last_object

    def get_list_member_with_smallest_index_greater_than_num(self, num, alist):
        for sem_scope_object in alist:
            if num <= sem_scope_object.struc_index:
                return sem_scope_object
        # in case alist is empty or num is greater than all of alist members' indices:
        return None

    def pprint(self):
        dur = 'Duration-like Strucs:\n    ' + '\n    '.join([sem_scope_object.pprint() for sem_scope_object in self.duration_like_strucs])
        tim = 'Timing-like Strucs:\n    ' + '\n    '.join([sem_scope_object.pprint() for sem_scope_object in self.timing_like_strucs])
        return dur + '\n' + tim


def create_sems(sig):
    """ For each parse, runs through parse.strucs to create an Instruction structure that represents the semantics of a maximal segment of strucs.

    AdminEvents represent atomic instructions, specifying single snapshot dose administration event. They don't have any information about repeating
    the event (e.g FREQ, Duration, Periodicity). If there is more than one AdminEvent in a Schedule, they have to have incompatible Timing properties (e.g. "take
    2 tabs in the am and 1 tab at bedtime" or be Specific_Day (which are all incompatible to each other).
    In the extreme, the non-initial AdminEvent may have only Timing info: "apply in the morning and also at bedtime".

    Schedules represent combinations of AdminEvents in time. They can specify several AdminEvents to take place during one day, and/or the freq of individual
    AdminEvent or group of AdminEvents, or their Periodicity, and/or their Duration (for how long should you continue those AdminEvents).
    If there are several Schedules in one Instruction, then they have to all either specify Duration or Calendar_Event that implies duration, such as "tomorrow",
    "now", "days 2-5", "day one").

    Instructions represent semantically maximal directions to patient on drug administration, one per drug. They may have several schedules, e.g. for tapering.


    We segment parse.strucs into Instructions/Schedules/AdminEvents using these semantic rules:
    - Rule 1 - Partitioning:
        - Instructions form a non-crossing partition of parse.strucs.
        - Schedules form a non-crossing partitition of parse.strucs of each Instruction of type Drug_Admin (except for those strucs that are proper semantic properties
          of Instruction, e.g. Indication -- these strucs are not assigned to any Schedule or AdminEvent)
        - AdminEvents form a non-crossing partition of all strucs of each Schedule that are properties of AdminEvent (not of Schedule or Instruction).

    - Rule 2 - Multiple SEMs:
        - If Instruction has more than 1 Schedule, each schedule has to have either Duration or Calendar_Event or TAPER.
        - If Schedule has more than 1 AdminEvent, each AdminEvent has to have Timing or Specific_Day.
        - Consecutive AdminEvents have to have incompatible Timings or Specific_Day (but see exceptions below). I.e. at least one of timing or Specific_Day strucs in AdminEvent_N
          has to be incompatible (can't coexist in the same AdminEvent) with at least one of the timing/Specific_Day strucs in AdminEvent_N+1.
        - Exceptions to the rule that each AdminEvent has to have an incompatible Timing or Specific_Day:
            -   Directive = "remove", e.g. "apply 1 patch for 12 hours then remove at bedtime". "REMOVE" should always start a new Event, because it represents a new
                action with a separate verb.
        - Exception to the rule that each Schedule has to have Duration or Calendar_Event or TAPER: If there is TAPER in a non-initial schedule, the left schedule doesn't need
          to have a Duration or Calendar_Event. E.g. "take one capsule at bedtime increase by 1 capsule every day". Taper (= "incrase by") has to start a new schedule.

    - Rule 3 - Boundaries:
        - In the case of multiple Schedules per one Instruction, we assume that the Duration struc is at the end (after all the strucs that properly
          belong to any of the AdminEvents that are part of this Schedule).
          Hence, we assume that the scope of the Duration struc is everything preceding it that does not already belong to a previous Schedule,
          as well as any subsequent strucs that are properties of schedule and that are compatible with the Schedule to which the Duration struc belongs to.

          The only exception we allow (in part to accomodate our dictionary) is to have the following form of Duration at the beginning:
          "Then for the next/first N days:", i.e. "Then_Chrono Duration:". In those cases, the scope of Schedule starts with this "right-facing" duration and
          ends with the following Then_Chrono or Duration, whichever comes first. (We therefore require that Then_chrono can't be part of any internal
          AdminEvents for such a Schedule. E.g., we outlaw "Then for the next 3 days: take 2 tabs in the AM, then 1 tab in the PM.").

        - The boundary between 2 consecutive Schedules starts with the first struc assignable to the later Schedule. The "assignable strucs" are those
          that:
            -- Are followed anywhere on the right by Duration
            and
            -- Are one of:
                --- Then_Chrono
                --- And_Conj
                --- any AdminEvent or Schedule struc that can't be assigned to the current event or schedule because it's slot is taken or is incompatible with that slot.
        - The boundary between 2 consecutive AdminEvents starts with the first struc that is already the boundary between Schedules or is assignable to the
          later Event. The "assignable strucs" are those that:
            -- Are followed anywhere on the right by Timing or Specific_Day
            and
            -- Are one of:
                --- Then_Chrono
                --- And_Conj
                --- any AdminEvent struc that can't be assigned to the current event because it's slot is taken or is incompatible with that slot.

        - If a struc is within a scope of a Schedule, then it can't start another Schedule. E.g., any Then_Chrono, And_Conj, etc. in the middle of a
          Schedule scope have to start an AdminEvent.

    PROBLEM1: sometimes the current Schedule doesn't have a DURATION -- it is implicitly or explicitly "from then on". E.g.,
    "take 3 tabs daily for 3 days, 2 tabs daily for 7 days, then take 1 tab daily."

    PROBLEM2: In our Dictionary, DURATION struc is at the start of most Schedules, not at the end. E.g., "Then for the next 3 days:"

    AHA: If you create a boundary between Schedules because of a presence of DURATION, then you can't later create other boundaries because of the same
    DURATION. So we create a function that maps coordiantes of Duration to the coordinates of the first struc compatible with that Duration, or None if we
    didn't encounter one yet.

    Algorithm:
    First, use Rule 2 to identify the centers of each Schedule and AdminEvent. Then create a map semcore2start that maps each Duration and Timing and Specific_Day and Calendar_Event
    struc to an index of its start in the parse.strucs array. This map is initialized by mapping all semcores to None, and then building up the actual values
    in the main sem loop.

    a simple heuristic: moving left to right:
    -   If the current struc is a sem property for some sem and can fill a semantic slot of the current sem we have been processing -- i.e. if that slot has
        not been filled (or if the slot can have multiple values, e.g. Timing or Freq, we check if the new struc is semantically compatible with existing fillers) --
        we simply assign it to that sem.

    -   Else: (i.e. the structure is a sem property for some sem but is incompatible with existing semantics, or it is a Special struc (Then_Chrono, And_Conj, Quant))
        we take it to signal the start of a new sem. Which one: AdminEvent, Schedule, or Instruction?
        - If the next Duration to the left is a right-scope Duration struc (e.g. "For the next 3 days:") then create a new AdminEvent.
        - Elif there is a Duration on the right that is unaccounted for by semcore2start, then we map it to this struc's index and create the start of a new Schedule.
          If the struc is of sem_type AdminEvent, also create a new AdminEvent.
        - Elif there is a Duration of the right that is accounted for by semcore2start (and the start of this Schedule is presumably on our left), then we start a
          new AdminEvent, not a new Schedule.
        - Elif (i.e. there is no Duration on the right) but the current Schedule has a Duration filled (ie. there is a Duration on the left) we assume that this
          is an implied Schedule to do something forever. E.g. "take 2 tabs daily for 10 days, then 1 tab daily." So we start a new Schedule with semantics "after that:".
        - Else: this is a start of a new AdminEvent.


        Problem: How do we account for:
        "Take 2 tabs in AM and 2 in PM for 3 days. Then for the next 7 days: take 1 tab in the AM _and_1_ in the PM. Then for the next 30 days: take 1 tab daily."
        How does the knowledge that the second Duration is right-scoped help us decide that "and 1 in the PM" is that start of a new AdminEvent in the second
        Schedule, not the start of a new Schedule.

    E.g. if the struc is a dose and the previous AdminEvent has a dose already, then this is a new AdminEvent from the existing Schedule.  -- Except
    that sometimes it's the start of a new AdminEvent and a new Schedule (e.g. "take 2 tabs for 2 days, 3 tabs for 3 days, and 4 tabs for 4 days".

    Complications arise when the struc does not correspond to any kind of semantic slot. These include Then_Chrono, And_Conj, Quant
    (which is likely to signal a clipped dose (e.g. "take 2 tabs in am, _3_ at noon, and 4 at bedtime").
    These strucs may signal the presence of either a new AdminEvent or a new Schedule. We decide which one using the following heuristic:
    If the struc is followed (not necessarilly immediately) by a Timing, then it is a new AdminEvent because a new
    AdminEvent has to differentiate itself from a previous one by some specification of timing. But if the struc is followed by Duration or by
    Calendar_Event (e.g. "tomorrow" or "Days 2-5") then it signals the presence of a new Schedule.

    In the case of strucs like Then_Chrono or And_Conj or Quant, even when they signal a new sem (AdminEvent or Schedule) they may not start it. E.g.
    "take 2 tabs now and then 1 tomorrow". Here, "then" and "1" indicate that a new Schedule has started, but only "and" starts that new schedule.  So how do
    we know if the new Schedule or AdminEvent is started or merely signaled by a struc of this type? If we the struc signals a new AdminEvent and
    the current_event has its Timing slot empty, then we don't need to start a new AdminEvent. If any of those slots have a value, we
    do start a new AdminEvent. Similarly, if the struc signals a new Schedule, then we check if the current_schedule has Duration slot and Calendar_Event
    empty. If they are both empty, then we don't need to start a new Schedule.

    """





    def create_new_sem(sem_type):
        """ Gets parameter sem_type in ('Instruction', 'Schedule', 'AdminEvent'). Returns updated tuple (current_inst, current_sched, current_evt)
        """

        current_inst = current_instruction
        current_sched = current_schedule
        current_evt = current_event

        if sem_type == 'Instruction':
            current_inst = parse.add_new_DrugAdmin()
            current_sched = None
            current_evt = None
        elif sem_type == 'Schedule':
            if not current_inst:
                # E.g., we are at the start of the sig.
                current_inst = parse.add_new_DrugAdmin()
            current_sched = current_inst.add_new_schedule()
            current_evt = None
        elif sem_type == 'AdminEvent':
            if not current_inst:
                current_inst = parse.add_new_DrugAdmin()
            if not current_sched:
                current_sched = current_inst.add_new_schedule()
            current_evt = current_sched.add_new_AdminEvent()

        return (current_inst, current_sched, current_evt)





    for parse in sig.parses:

        sem_scope_data = SemScopeData(parse.strucs)
        #print sem_scope_data.pprint()

        current_instruction = None
        current_schedule = None
        current_event = None

        for (struc_index, struc) in enumerate(parse.strucs):
            duration_on_the_left = sem_scope_data.get_list_member_with_largest_index_less_than_num(struc_index, sem_scope_data.duration_like_strucs)
            duration_on_the_right = sem_scope_data.get_list_member_with_smallest_index_greater_than_num(struc_index, sem_scope_data.duration_like_strucs)
            timing_on_the_right = sem_scope_data.get_list_member_with_smallest_index_greater_than_num(struc_index, sem_scope_data.timing_like_strucs)

            if struc.label == 'THEN_CHRONO':
                # "Then" signals the start of a new Schedule or a new AdminEvent.
                # It signals new Schedule if it is followed by a new duration (or by most types of calendar_event). So we use a lookahead loop to verify this.
                # E.g., "2 tabs now then 1 tab every morning for the next 4 days", or " then 1 tab tomorrow" or "then 1 tab on Monday"
                # However, "then" can signal merely a new AdminEvent if it consists in change of timing from the previous AdminEvent of the same schedule.
                # E.g. �take 2 tabs in the morning then 1 in the evening then 1 at bedtime for pain as needed�


                # To determine if we need to start a new Schedule:
                # Check if there is Duration to the right with unassigned start.
                #       If yes, create a new Schedule and assign this struc as the start of the schedule in schedule_start_dict.
                #       If schedule_start_dict has a positive index assigned to the Schedule, do not create a new schedule. Proceed to check if it is the start of a new AdminEvent
                #       If there are no Duration strucs to the right or if the next duration_on_the_right has None as it's left_scope_start, then:
                #           a. Check if we are at the start of a left-scoped Duration, i.e. if we are at "Then, Duration:". If so, start a new Schedule.
                #           b. Else check if there is a Duration to the right with left_scope_start assigned to it. If so, do not create a new schedule.
                #              Proceed to check if it is the start of a new AdminEvent
                #           c. If there is no Duration on the right but the current Schedule has a Duration filled we assume that this
                #               is an implied Schedule to do something forever. E.g. "take 2 tabs daily for 10 days, then 1 tab daily." So we start a new Schedule with semantics "after that:".

                if duration_on_the_right and duration_on_the_right.is_right_scoped and duration_on_the_right.struc_index == struc_index + 2:
                    # I.e. we are the start of "Then for the next n days: ..."
                    (current_instruction, current_schedule, current_event) = create_new_sem('Schedule')
                    struc.accounted_for_by_sem = current_schedule
                elif duration_on_the_right and duration_on_the_right.is_left_scoped and duration_on_the_right.left_scope_start is None:
                    duration_on_the_right.left_scope_start = struc_index
                    (current_instruction, current_schedule, current_event) = create_new_sem('Schedule')
                    struc.accounted_for_by_sem = current_schedule
                elif duration_on_the_right is None and duration_on_the_left:
                    # If there is no Duration on the right but the current Schedule has a Duration filled we assume that this
                    # is an implied Schedule to do something forever. E.g. "take 2 tabs daily for 10 days, then 1 tab daily." So we start a new Schedule with semantics "after that:".
                    (current_instruction, current_schedule, current_event) = create_new_sem('Schedule')
                    struc.accounted_for_by_sem = current_schedule
                elif current_schedule is None:
                    # A new AdminEvent can only be started by a Then_Chrono if it is in the middle of a Schedule, i.e. a prior AdminEvent has well been started.
                    # So we are probably at the beginning of a new Instruction. Start a new Schedule then.
                    (current_instruction, current_schedule, current_event) = create_new_sem('Schedule')
                    struc.accounted_for_by_sem = current_schedule
                elif timing_on_the_right:
                    # We didn't start a new Schedule, so Then_Chrono has to signify a start of a new Admin_Event
                    # Since we are starting a new AdminEvent, the scope of Timing on the right has to start here.
                    timing_on_the_right.left_scope_start = struc_index
                    (current_instruction, current_schedule, current_event) = create_new_sem('AdminEvent')
                    struc.accounted_for_by_sem = current_event
                elif current_event is None:
                    (current_instruction, current_schedule, current_event) = create_new_sem('AdminEvent')
                    struc.accounted_for_by_sem = current_event

            elif struc.label == 'AND_CONJ':
                # "AND" could signal the start of a new Schedule or a new AdminEvent.
                # It signals new Schedule if it is followed by a new duration or calendar_event. So we use a lookahead loop to verify this.
                # E.g., "2 tabs day 1 and 3 tabs days 2-5", or "2 tabs today and 1 tab tomorrow" or "2 tabs now and 1 tab on Monday"
                # However, "and" can signal merely a new AdminEvent if it consists in change of timing from the previous AdminEvent of the same schedule.
                # E.g. �take 2 tabs in the morning and 1 in the evening�

                # EXCEPT you don't need incompatible Timing to start a new AdminEvent! All you need is "AND_Chrono" before Timing.
                # E.g. "Take 2 tabs every 3 hours and at bedtime." AND_Chrono there starts a new AdminEvent.

                if duration_on_the_right and duration_on_the_right.is_left_scoped and duration_on_the_right.left_scope_start is None:
                    duration_on_the_right.left_scope_start = struc_index
                    (current_instruction, current_schedule, current_event) = create_new_sem('Schedule')
                    struc.accounted_for_by_sem = current_schedule
                elif current_schedule is None:
                    # A new AdminEvent can only be started by a Then_Chrono if it is in the middle of a Schedule, i.e. a prior AdminEvent has well been started.
                    (current_instruction, current_schedule, current_event) = create_new_sem('Schedule')
                    struc.accounted_for_by_sem = current_schedule
                elif timing_on_the_right:
                    # We are starting a new AdminEvent, so break the scope of timing_on_the_right to here.
                    timing_on_the_right.left_scope_start = struc_index
                    (current_instruction, current_schedule, current_event) = create_new_sem('AdminEvent')
                    struc.accounted_for_by_sem = current_event
                elif current_event is None:
                    (current_instruction, current_schedule, current_event) = create_new_sem('AdminEvent')
                    struc.accounted_for_by_sem = current_event

            elif DrugAdmin.is_valid_property(struc.label):
                if current_instruction and not struc.is_semantically_incompatible_with_given_sem(current_instruction):
                    # If the struc is compatible with an existing current_instruction, just add this struc to that current_instruction.
                    current_instruction.add_property_value(struc.label, struc)
                elif not current_instruction:
                    # Need to create the first Instruction and assign the struc to it.
                    (current_instruction, current_schedule, current_event) = create_new_sem('Instruction')
                    current_instruction.add_property_value(struc.label, struc)
                #else:
                # We don't create a new instruction because of incompatibility, because all the cases we have seen (e.g. multiple AS_NEEDED) are just stupid repetitions
                # which don't intend to signal new instruction.


            elif Schedule.is_valid_property(struc.label):
                if duration_on_the_left and duration_on_the_left.is_right_scoped and duration_on_the_left.struc_index == struc_index:
                    # I.e. we are the start of a right-scoped Duration: "For the next n days: ...". Need to start a new Schedule.
                    (current_instruction, current_schedule, current_event) = create_new_sem('Schedule')
                elif duration_on_the_left and duration_on_the_left.is_right_scoped and duration_on_the_left.right_scope_end >= struc_index:
                    # We are in the scope of an already-started right-scoped Duration Schedule, so don't start anything new.
                    pass
                elif duration_on_the_right and duration_on_the_right.is_left_scoped and duration_on_the_right.left_scope_start is None:
                    duration_on_the_right.left_scope_start = struc_index
                    (current_instruction, current_schedule, current_event) = create_new_sem('Schedule')
                elif current_schedule is None:
                    (current_instruction, current_schedule, current_event) = create_new_sem('Schedule')
                elif struc.label == 'TAPER':
                    # Taper always starts a new schedule. E.g. "take one capsule at bedtime increase by 1 capsule every day"
                    duration_on_the_right.left_scope_start = struc_index
                    (current_instruction, current_schedule, current_event) = create_new_sem('Schedule')
                current_schedule.add_property_value(struc.label, struc)

            elif AdminEvent.is_valid_property(struc.label):
                if current_event and not struc.is_semantically_incompatible_with_given_sem(current_event):
                    # If the struc is compatible with an existing current_event, just add this struc to that current_event.
                    current_event.add_property_value(struc.label, struc)
                else:
                    start_new_event = True
                    if duration_on_the_left and duration_on_the_right and duration_on_the_right.is_left_scoped:
                            # This may be a multi_schedule instruction which doesn't have Timing in the first schedule's events.
                            # E.g. "take 7 tab twice daily for 5 days, 2 in am & 3 in pm for 4 days,then one tab once daily"
                            start_new_event = True
                    elif current_event and not current_event.timing:
                        # We are in a bind. We have an incompatible structure to the current_event, so on the one hand we need to close that event and
                        # start a new one. But on the other hand, current_event has not Timing, so it is not really an Event.
                        # The situation is typically due to repeated directive, e.g. "use 1 apply vaginally every night at bedtime"
                        start_new_event = False
                        if struc.label == 'DIRECTIVE' and current_event.directive:
                            if struc.value == 'remove':
                                # "Remove" should always start a separate Event for it represents a separate action with
                                # it's own properties. The previous event doesn't have to have Timing. E.g. "apply 1 patch for 10 hours, then remove at bedtime"
                                # For the frequent case of Duration on/off (e.g. "apply 1 patch every day on 12 hours off for 12 hours") we DO want to create a new
                                # Schedule (and new Event to start it) because we have 2 Durations to deal with. We finagle the problem of this being really an intra-day duration in
                                # process_special_duration_cases()
                                # e.g. "apply one patch once weekly as directed for three weeks. leave off for one week then repeat cycle". "remove" with Duration should start a new event.
                                start_new_event = True
                            elif (
                                    (current_event.directive.value == 'use' and struc.value in ('inject', 'instill', 'mix', 'dissolve', 'chew', 'inject', 'take')) or
                                    (current_event.directive.value == struc.value) or
                                    (current_event.directive.value == 'use' and struc.value == 'rinse' and (not (current_event.dose) or current_event.dose.form.value in ('capful', 'ounce', 'ml', 'teaspoon', 'tablespoon')))
                                ):
                                # just make this directive be the event directive.
                                # The case of rinse is this: we want to cover cases such as "use 1 capful to rinse.." but we want to avoid "use 1 puff twice a day, rinse mouth after"
                                current_event.directive.rules_used.append('*removed_from_sem_in_create_new_sem*')
                                current_event.remove_property('directive', current_event.directive)
                                current_event.add_property_value(struc.label, struc)
                                start_new_event = False
                            elif current_event.directive.value in ('inject', 'instill', 'mix', 'dissolve', 'chew', 'inject', 'take') and struc.value in ('use', 'take'):
                                # remove this directive. E.g. "chew one tablet take one tablet three times a day before meals "
                                struc.rules_used.append('*removed_from_sem_in_create_new_sem*')
                                start_new_event = False
                            elif current_event.directive.value in ('mix', 'dissolve') and struc.value in ('drink', 'inject'):
                                # This is a multi-event situation, e.g. "mix with 1 cc diluent and inject intramuscularly"
                                # We will remove "ALSO:" from the transduction in process_schedule().
                                start_new_event = True
                            elif struc.value in ('stop'):
                                # This is a multi-schedule situation, e.g. "take 1 tablet by mouth once daily for 21 days, skip 7 days and repeat cycle"
                                start_new_event = True
                        elif struc.label == 'DOSE' and current_event.dose:
                            if current_event.dose.form.value != struc.form.value:
                                # E.g. "take 1 tablet and 1 gelcap by mouth daily". So even though there is no timing difference, we can pretend that we first take the tablet then take the capsule.
                                start_new_event = True
                            elif current_event.specific_day:
                                # we should start a new event because this is likely to be a switch of doses b/c of change in days:
                                # eg. "take one tablet by mouth on monday thru friday one and half tablet saturday and sunday"
                                start_new_event = True
                            else:
                                # probably repeat "take 1 tablet 1 tablet", so ignore the second Dose
                                start_new_event = False
                        elif struc.label == 'SPECIFIC_DAY' and current_event.specific_day:
                            start_new_event = True

                    if start_new_event:
                        # Need to create a new event and assign the struc to it.
                        # But first check if you need to start a new Schedule.
                        if duration_on_the_right and duration_on_the_right.is_left_scoped and duration_on_the_right.left_scope_start is None:
                            duration_on_the_right.left_scope_start = struc_index
                            (current_instruction, current_schedule, current_event) = create_new_sem('Schedule')
                        elif current_schedule is None:
                            (current_instruction, current_schedule, current_event) = create_new_sem('Schedule')
                        (current_instruction, current_schedule, current_event) = create_new_sem('AdminEvent')
                        if timing_on_the_right:
                            # We are starting a new AdminEvent, so the scope of Timing on the right has to start here.
                            timing_on_the_right.left_scope_start = struc_index
                        current_event.add_property_value(struc.label, struc)


        parse.sem_scope_data = sem_scope_data

        #print sem_scope_data.pprint()
        #s = parse.show_struc_assignment_to_sems(include_struc_details = True, omit_coords = False)
        #print s


def modify_sems(sig):
    """  Makes appropriate deductions and modifications to Sems.
    """

    def modify_directive_and_form(instruction, directive, form, route, site, substrate):

        if site and site.value in ('eye', 'ear'):
            if directive and directive.value in ('instill', 'apply', 'place', 'use') and not form:
                form = Form(form_name = 'drop', plurality = '', constituents = [])
                form.rules_used.append('*deduced_in_modify_directive_and_form_01*')
            elif not directive and not form:
                form = Form(form_name = 'drop', plurality = '', constituents = [])
                form.rules_used.append('*deduced_in_modify_directive_and_form_02*')
            elif form and form.value == 'application' and (not directive or directive.value in ('instill', 'apply', 'place', 'use')):
                form = Form(form_name = 'drop', plurality = '', constituents = [])
                form.rules_used.append('*deduced_in_modify_directive_and_form_03*')
        elif site and site.value == 'nostril':
            if not form and (not directive or directive.value in ('instill', 'apply', 'use')):
                form = Form(form_name = 'spray', plurality = '', constituents = [])
                form.rules_used.append('*deduced_in_modify_directive_and_form_04*')
                if directive and directive.value != 'instill':
                    directive.value = 'spray'
                    directive.rules_used.append('*deduced_in_modify_directive_and_form_05*')
            elif form and form.value == 'drop':
                if directive and directive.value != 'instill':
                    directive.value = 'instill'
                    directive.rules_used.append('*deduced_in_modify_directive_and_form_06*')
                elif not directive:
                    directive = Directive(value = 'instill')
                    directive.rules_used.append('*deduced_in_modify_directive_and_form_07*')


        if route:
            if not form and directive and route.value == 'orally':
                # oral route maps to form = tablet. But avoid it for things like 'rinse'
                if directive.value not in ('rinse', 'spray', 'mix', 'dissolve'):
                    form_name = Route.value_to_likely_form.get(route.value, '')
                    if form_name:
                        form = Form(form_name = form_name, plurality = '', constituents = [])
                        form.rules_used.append('*deduced_in_modify_directive_and_form_08*')
            elif not form:
                # e.g. "use 1 vaginally every night"
                form_name = Route.value_to_likely_form.get(route.value, '')
                if form_name:
                    form = Form(form_name = form_name, plurality = '', constituents = [])
                    form.rules_used.append('*deduced_in_modify_directive_and_form_09*')
            elif 'deduced' in ''.join(form.rules_used):
                form_name = Route.value_to_likely_form.get(route.value, '')
                if form_name:
                    form.value = form_name
                    form.rules_used.append('*deduced_in_modify_directive_and_form_10*')

            if not directive and form and form.value in ('spray', 'puff'):
                # E.g. "1 puff by mouth twice a day". Oral should not imply take but "use" or "inhale"
                if route.value == 'orally':
                    directive = Directive('inhale')
                else:
                    directive = Directive('use')
                directive.rules_used.append('*deduced_in_modify_directive_and_form_11*')
            elif not directive:
                # e.g. "1 ORAL EACH DAY" -> then directive is "take", or "1 topical 3 times a day" -> directive is 'apply'
                directive_value = Route.value_to_likely_directive.get(route.value, '')
                if directive_value:
                    directive = Directive(directive_value)
                    directive.rules_used.append('*deduced_in_modify_directive_and_form_07*')
            elif directive.value == 'use' and route.value in ('orally', 'vaginally') and form.value not in ('spray', 'puff', 'vial'):
                # avoid "use 1 spray by mouth", because here "use" is not "take"
                directive_value = Route.value_to_likely_directive.get(route.value, '')
                if directive_value:
                    if directive.accounted_for_by_sem and directive.accounted_for_by_sem.anaphora:
                        # don't modify the verb in the dictionary. This is a hack to let "Use this medicine by mouth" be processed as "use", not as "take"
                        pass
                    else:
                        directive.value = directive_value
                        directive.rules_used.append('*deduced_in_modify_directive_and_form_08*')

        if substrate and directive and directive.value == 'take' and form and form.value in ('capful', 'gram', 'packet', 'scoop'):
            # "take 1 capfull in 8 oz of water or juice and drink" means "mix"
            # But "take 1 tablet with a glass of water" or (from dictionary: "Take this medicine with water.") should not be changed to "mix"
            directive.value = 'mix'
            directive.rules_used.append('*deduced_in_modify_directive_and_form_09*')

        if not directive and not form:
            # either there is no route or we can't figure out what the directive should be from the route
            # e.g. "ONCE DAILY"
            directive = Directive('use')
            directive.rules_used.append('*deduced_by_defaulting_to_use_bc_of_absence_of_form_or_route*')

        if form and not directive:
            directive_name = Form.value_to_likely_directive.get(form.value, '')
            if directive_name:
                directive = Directive(directive_name)
                directive.rules_used.append('*deduced_in_modify_directive_and_form_10*')
        if directive and not form:
            form_name = Directive.directive_to_likely_form.get(directive.value, '')
            if form_name:
                form = Form(form_name = form_name, plurality = '', constituents = [])
                form.rules_used.append('*deduced_in_modify_directive_and_form_11*')


        if directive and form:
            if form.value == 'drop' and directive.value in ('place', 'apply') and site and site.value != 'knees':
                # Use "instill N drops" instead of "place N drops"
                # But not in the case of site = knees, because there is an expression "apply 40 drops to both knees" that
                # is used for PENNSAID, an NSAID topical solution for osteoarthritis of the knee.
                directive.value = 'instill'
                directive.rules_used.append('deduced_in_modify_directive_and_form_12*')
            elif form.value == 'ring' and directive.value == 'place':
                # Use "insert ring vaginally" instead of "place ring"
                directive.value = 'insert'
                directive.rules_used.append('*deduced_in_modify_directive_and_form_13*')
            elif form.value == 'patch' and directive.value == 'use':
                # "apply 1 patch" instead of "use 1 patch"
                directive.value = 'apply'
                directive.rules_used.append('*deduced_in_modify_directive_and_form_14*')
            elif form.value == 'inhalation' and directive.value == 'use':
                directive.value = 'inhale'
                directive.rules_used.append('*deduced_in_modify_directive_and_form_15*')
            elif form.value == 'ounce' and directive.value in ('take', 'give'):
                directive.value = 'drink'
                directive.rules_used.append('*deduced_in_modify_directive_and_form_16*')
            elif form.value == 'ounce' and directive.value in ('swish', 'swish and spit'):
                directive.value = 'rinse'
                directive.rules_used.append('*deduced_in_modify_directive_and_form_16*')
            #elif form.value == 'puff' and directive.value == 'take':
            #    # "inhale 1 puff" instead of "take 1 puff"
            #    directive.value = 'inhale'


        return (directive, form)


    def modify_directive_from_route(directive, route, form):
        """ Used twice: for each event and also for instruction as a whole to determine the primary_directive
            because sometimes there is no directive in the event but it can be deduced (in determine_primary_directive_and_form())
        """

        if directive and route:
            expected_directive_value = Route.value_to_likely_directive.get(route.value, '')
            if expected_directive_value and expected_directive_value != directive.value:
                # "DISSOLVE 1 TABLET SUBLINGUALLY " instead of "USE 1 TABLET SUBLINGUALLY "
                # or "apply topically" instead of "use topically"

                # but we are being cautious now because this is changing the event itself. For now, just go with "sublingually"
                if route.value == 'sublingually' and (not form or form.value != 'spray'):
                    # Avoid "use 1 spray under the tongue as needed" because we don't dissolve sprays.
                    directive.value = expected_directive_value
                    directive.rules_used.append('*deduced_in_modify_directive_from_route*')
                elif route.value == 'intramuscularly':
                    directive.value = expected_directive_value
                    directive.rules_used.append('*deduced_in_modify_directive_from_route*')



    def determine_primary_directive_and_form(instruction):
        form = None
        directive = None
        directive_from_as_directed = None  # Directive deduced from "AS_DIRECTED" phrases e.g. "instill as directed by doctor"
        route = None
        site = None
        substrate = None
        for schedule in instruction.schedules:
            for event in schedule.events:
                if not directive and event.directive and event.directive.value not in ('stop'):
                    # "STOP" can never be a primary directive. You have to stop something first.
                    directive = event.directive
                elif not directive and event.dir_attribute and 'per_sliding_scale' in ''.join([attr.value for attr in event.dir_attribute]):
                    directive = Directive('inject')
                    directive.rules_used.append('*deduced_from_dir_attribute_sliding_scale*')
                elif directive and event.directive:
                    if directive.value == 'mix' and event.directive.value == 'inject':
                        # When there are several directives (from different events), sometimes the first one is not the main one. E.g. "mix with 1 cc diluent and inject im..."
                        directive = event.directive
                    elif directive.value == 'use' and event.directive.value in ('apply', 'take', 'instill', 'insert'):
                        # the first directive, "use" is just a place holder. E.g. "use 1 apply vaginally every night at bed time". Replace it with the more informative one.
                        directive = event.directive
                if not form and event.dose and event.dose.form.value:
                    form = event.dose.form
                if not route and event.route:
                    route = event.route
                if not site and event.site:
                    site = event.site
                if not substrate and event.substrate:
                    substrate = event.substrate
            # If there is no directive in any Event, see if you can use the implied directive in AS_DIRECTED (e.g. "inject as directed twice a day"
            if (not directive and not directive_from_as_directed and schedule.as_directed and
                schedule.as_directed.directive and schedule.as_directed.directive in Directive.permissible_values):
                directive_value = schedule.as_directed.directive
                directive = Directive(directive_value)
                directive.rules_used.append('*deduced_determine_primary_directive_and_form_0*')
        if not directive and directive_from_as_directed:
            directive = directive_from_as_directed

        if not form:
            # sometimes the FORM is not parse into a dose, and therefore is not part of the Sem hierarchy.
            parse = instruction.parse
            form_strucs = [struc for struc in parse.strucs if struc.label == 'FORM']
            if form_strucs:
                form = form_strucs[0]

        # Modify directive and form
        (directive, form) = modify_directive_and_form(instruction, directive, form, route, site, substrate)
        modify_directive_from_route(directive, route, form)

        if form and directive and not route:
            directive_is_plausible_for_oral_route = directive.value in ('take', 'give', 'chew', 'chew and swallow', 'dissolve', 'inhale')
            form_is_plausible_for_oral_route = form.value in ('tablet', 'capsule', 'puff', 'lozenge', 'ml', 'teaspoon', 'tablespoon', 'dropperful', 'puff', 'unit', 'packet')
            if directive_is_plausible_for_oral_route and form_is_plausible_for_oral_route:
                route = Route(value = 'orally', constituents = [])
                route.rules_used.append('*deduced_determine_primary_directive_and_form_6*')
        instruction.primary_directive = directive
        instruction.primary_form = form
        instruction.primary_route = route

    def insert_directive_into_event(event):
        """ Usefull for processing very typical things like "1 TABLET ONCE DAILY"        """

        if not event.directive and event.schedule.instruction.primary_directive and (event.dose or event.route or event.site):
            event.directive = event.schedule.instruction.primary_directive
            if not event.directive.accounted_for_by_sem:
                event.directive.accounted_for_by_sem = event

    def insert_substrate_into_event(event):
        """ Usefull for processing very typical things like "mix and drink by mouth once a day"
            We insert substrate "mix in water" if it is missing.
            But need to be careful to avoid dissolving lozenges under the tongue and dissolving stuff in water for the purpose of rinsing/swishing
        """


        if event.directive and event.directive.value in ('mix', 'dissolve') and not event.substrate:
            # Avoid "dissolve 1 tab under the tongue"
            if (
                    (event.schedule.instruction.primary_form and event.schedule.instruction.primary_form in ('tablet', 'lozenge')) or
                    (event.dose and event.dose.form.value in ('tablet', 'lozenge', 'drop')) or      # reason for excluding drop: "'mix 4 drops in 4 ounces of water . swish and spite three times daily'"
                    (event.route and event.route.value == 'sublingually')
                ):
                return None
            substrate = Substrate(value = 'water')
            substrate.rules_used.append('*deduced_in_insert_substrate_into_event*')
            event.substrate = substrate


    def modify_directive_in_event(event):
        """ Modifies directive in several cases:
            -   When the source says "off for 1 week" we transduce it as directive "REMOVE"
                because most of the time it is about patches or rings that need to be removed for a 1 week or for 12 hours.
                But for other things such as tablets, applicatorfulls, etc. 'off' really means stop.
            -   Dissolve packets in water or take 17 grams with 8 oz of liquid need to have direcitve = mix.
        """

        if event.directive and event.directive.value == 'remove' and 'remove' not in ''.join(event.directive.constituents):
            form = event.schedule.instruction.primary_form
            if form and form.value not in ('ring', 'patch'):
                event.directive.value = 'stop'
                event.directive.rules_used.append('*deduced_in_modify_directive_in_event*')
        elif event.substrate and event.dose and event.dose.form.value in ('packet', 'ounce', 'gram', 'ml', 'liter', 'cc', 'scoop'):
            if event.directive and event.directive.value != 'mix':
                event.directive.value = 'mix'
                event.directive.rules_used.append('*deduced_from_substrate*')
            else:
                directive = Directive('mix')
                directive.rules_used.append('*deduced_from_substrate*')
                event.add_property_value('directive', directive)
            event.schedule.instruction.primary_directive = event.directive


    def insert_dose_into_event_from_quant(parse):
        """ Missing dose: If QUANT is in the parse.strucs (i.e. it hasn't been fully parsed and there is no dose in the event to which
            QUANT presumably (positionally) belongs, assume that QUANT is part of DOSE and insert DOSE.

            Typical examples:
                1 ORAL EACH DAY
                take 1 BID
                1TID
                1 bid for 7 days then 2 tid

            But avoid "QUANT times", e.g. "inject intramuscularly once"

            The Event that Quant belongs it: if the event closely follows the Quant, then use that event; otherwise use the event
            that last precedes Quant.
        """

        def find_quant_and_event(starting_struc_index):
            last_event = None
            last_schedule = None
            quant = None
            quant_position = None

            # First, find the Quant, if any, and the preceding Event, if any (to which the Quant might belong as part of Dose).
            for struc_num, struc in enumerate(parse.strucs[starting_struc_index:]):
                sem = struc.accounted_for_by_sem
                if sem and isinstance(sem, AdminEvent):
                    last_event = sem
                    last_schedule = last_event.schedule
                elif sem and isinstance(sem, Schedule):
                    last_schedule = sem
                elif struc.label == 'QUANT':
                    quant = struc
                    quant_position = struc_num
                    break

            if quant is None:
                return None

            # make sure that we are not seeing "QUANT times"
            tail_string = get_struc_labels(parse.strucs[quant_position:], delimiter = '', omit_spaces_and_punctuation = False)
            times_pattern = re.compile('^QUANT\s*time')
            if times_pattern.search(tail_string):
                return None


            # Now find the first event that follows Quant in the next 3 positions (we don't want to
            # search further than that because it is likely then that it is the wrong Event. E.g. in "Take 1 bid first day, days 1-7 take 2 tid"
            # the first "1" clearly belongs to the first Event, even though the last struc of the first Event is the first word of the sig (first "take").
            # But in "1 orally bid" Quant belongs to the succeeding event.
            succeeding_event = None
            search_through_position = (starting_struc_index + quant_position + 4) if last_event else (len(parse.strucs) + 1)
            for struc in parse.strucs[starting_struc_index + quant_position + 1:search_through_position]:
                sem = struc.accounted_for_by_sem
                if struc.label in ('THEN_CHRONO', 'AND_CONJ', 'DOSE', 'FORM', 'QUANT'):
                    # We are definitely in a new Event to which the Quant doesn't belong. It belongs to a previous Event, if any.
                    break
                elif sem and not last_schedule and isinstance(sem, Schedule):
                    last_schedule = sem
                elif sem and isinstance(sem, AdminEvent):
                    succeeding_event = sem
                    if not last_schedule:
                        last_schedule = succeeding_event.schedule
                    break

            if succeeding_event:
                event_to_insert_in = succeeding_event
            elif last_event:
                event_to_insert_in = last_event
            elif last_schedule and not last_schedule.events:
                # We have no events before the quant and the stuff to the right is either not in any event (e.g. "1 bid" is just QUANT FREQ, hence no Events in the Sem)
                # or it is a different event (e.g. "1 bid for 7 days then 1 tid". So we need to create a new event in the current schedule.
                event_to_insert_in = last_schedule.add_new_AdminEvent()
            else: # There are no events before or after Quant, and no schedule. We just have some lonely number here. Quit.
                return None

            if event_to_insert_in.dose:
                return None

            instruction = event_to_insert_in.schedule.instruction
            #determine_primary_directive_and_form(instruction)
            dose_plurality = 'singular' if quant.value == 1 else 'plural'

            # There is no explicit directive if the event has no directive and the directive for the instruction has been deduced, not explicitly given in some other event
            no_explicit_directive = (not event_to_insert_in.directive and
                                     (not instruction.primary_directive or 'deduced_by_defaulting_to_use_bc_of_absence_of_form_or_route' in ''.join(instruction.primary_directive.rules_used)))
            no_explicit_form = (not instruction.primary_form or 'deduced' in ''.join(instruction.primary_form.rules_used))

            if no_explicit_directive and no_explicit_form:
                directive = Directive('take')
                directive.rules_used.append('*deduced_by_defaulting_to_take_in_presence_of_only_quantifier')
                instruction.primary_directive = directive
                form = Form(form_name = 'tablet', plurality = dose_plurality, constituents = [])
                form.rules_used.append('*deduced_by_defaulting_to_tablet_in_presence_of_only_quantifier*')
                if not instruction.primary_form:
                    # plurality can change if there are multiple events in the same instruction, e.g. "1bid for 7 days then 2 tid",
                    # so don't change the plurality of primary_form when you come to the second event
                    instruction.primary_form = form
            elif not instruction.primary_form:
                # That means primary_directive exists or was deduced already
                directive = event_to_insert_in.directive if event_to_insert_in.directive else instruction.primary_directive
                form_name = Directive.directive_to_likely_form.get(directive.value, '')
                if not form_name:
                    return None
                form = Form(form_name, plurality = dose_plurality, constituents = [])
                form.rules_used.append('*deduced_from_directive_in_presence_of_only_quantifier*')
                instruction.primary_form = form
            elif instruction.primary_form.plurality == dose_plurality:
                form = instruction.primary_form
            elif not instruction.primary_form.plurality or instruction.primary_form.plurality == 'plurality_either':
                instruction.primary_form.plurality = dose_plurality
                form = instruction.primary_form
            else:
                form = Form(form_name = instruction.primary_form.value, plurality = dose_plurality, constituents = [])
                form.rules_used.append('*deduced_from_primary_form_in_presence_of_only_quantifier*')

            if form:
                dose = Dose(quant, form)
                event_to_insert_in.add_property_value('DOSE', dose)
                quant.rules_used.append('*assigned_to_dose_in_insert_dose_into_event_from_quant*')
                return quant_position

        number_of_quants = parse.flatten().count('QUANT')
        starting_struc_index = 0
        for quant_occurence_num in range(number_of_quants):
            quant_position = find_quant_and_event(starting_struc_index = starting_struc_index)
            if not type(quant_position) is int:
                break
            else:
                starting_struc_index = quant_position + 1



    def modify_nightly_bedtime_timings(event):
        """ Remove redundant Timing(landmark = 'night') if it is conjoined with Timing(landmark = 'beditme')

            E.g. "take nightly at bedtime" should have freq=every day and only 1 timing ('bedtime')
        """

        if event.timing and len(event.timing) > 1:
            night_struc = None
            bedtime_struc = None
            for struc in event.timing:
                if struc.landmark == 'night':
                    night_struc = struc
                if struc.landmark == 'bedtime':
                    bedtime_struc = struc
            if bedtime_struc and night_struc:
                event.remove_property('timing', night_struc)
                pos = parse.strucs.index(night_struc)
                if len(parse.strucs) > pos + 1 and parse.strucs[pos + 1].label == ' ':
                    # remove the space after Timing
                    parse.strucs[pos:pos+2] = []
                else:
                    parse.strucs.remove(night_struc)

    def modify_specific_day(event):
        """ Remove "every day" in Specific_Day struc, because it is a contradiction "take 1 tab daily Mon-Fri" doesn't mean "take every day", only every day in that range.

        """

        if event.specific_day:
            schedule = event.schedule
            remove_strucs = []
            if schedule.periodicity:
                for periodicity in schedule.periodicity:
                    if periodicity.time_unit.value == 'day':
                        remove_strucs.append(periodicity)
            for struc in remove_strucs:
                    schedule.periodicity.remove(struc)
                    struc.accounted_for_by_sem = None
                    struc.rules_used.append('*removed_from_sem_in_modify_specific_day*')

    def remove_and_conj_between_indications(parse):
        """ Removes the AND_CONJ struc that is unaccounted for by sems if it is between 2 indications because in that case it doesn't indicate new AdminEvent or new Schedule onset.
        """

        remove_strucs = []
        for i, struc in enumerate(parse.strucs):
            if struc.label == 'AND_CONJ' and not struc.accounted_for_by_sem:
                if len(parse.strucs) > i + 2 and parse.strucs[i+1].is_space_or_punctuation_only() and parse.strucs[i+2].label == 'INDICATION':
                    remove_strucs += parse.strucs[i:i+2]
        parse.strucs = [struc for struc in parse.strucs if struc not in remove_strucs]

    def remove_strucs_following_unparsed_maxdose(parse):
        """ If there is a maxdose that is unparsed (there is a loose "max" string still remaining in the strucs, then remove the subsequent strucs from the sem
            because they create an illusion of multiple Events.
        """

        remove_after = None
        for i, struc in enumerate(parse.strucs):
            if 'max' in struc.label:
                remove_after = i
            if remove_after is not None and struc.accounted_for_by_sem:
                sem = struc.accounted_for_by_sem
                prop = struc.label.lower()
                if prop in sem.__dict__:
                    sem_prop_value = sem.__dict__[prop]
                    if type(sem_prop_value) == list:
                        sem_prop_value.remove(struc)
                    else:
                        sem.__dict__[prop] = None

    # NOT USED. Consider REMOVING 04/04/2013. Problem: Sometimes there are no events in the last schedule, just "after that: take daily or every other day". E.g. "1bid for 7 days then 2 tid"
    def remove_extra_schedules(instruction):
        """ If there is a non-initial schedule that has no events, remove it so as not to create extra "After that:"
            One exception is the special "on_off" durations, e.g. Apply 12 hours on 12 hours off. We treat them as durations syntactically but then
            remove them during atom match processing (process_special_duration_cases())
        """

        #if instruction.parse.is_fully_parsed(): return

        bad_schedules = []
        for sched_num, schedule in enumerate(instruction.schedules):
            if not schedule.events and sched_num > 0 and not (schedule.duration and schedule.duration.on_off):
                bad_schedules.append(schedule)

        if bad_schedules:
            for struc in instruction.parse.strucs:
                if struc.accounted_for_by_sem in bad_schedules:
                    struc.accounted_for_by_sem = None
        for schedule in bad_schedules:
            instruction.schedules.remove(schedule)

    def remove_duplicate_indications(instruction):
        """ Example: "take 150 ml by mouth twice daily for constipation for bowel movement"
        """

        bad_indications = []
        conditions = set()
        for indication in instruction.indication:
            if indication.condition in conditions:
                bad_indications.append(indication)
            else:
                conditions.add(indication.condition)
        for indication in bad_indications:
            instruction.remove_property('indication', indication)
            indication.rules_used.append('*removed_from_sem_in_remove_duplicate_indications*')

    def remove_specious_events(sched_num, schedule):
        """ If there is an event that has only 1 property that makes no sense on it's own, it is probably a mistake in the source or in parse.
            We don't want to create specious "Also:", so remove the event

            Example:
                SITE: use 2 sprays in each nostril twice a day as directed in each nostril
                ROUTE: "take one tablet by mouth daily and by mouth 2 times a day"
                STOP directive: "take one tablet daily stop tekturna" or "take one tablet daily stop 600mg". The legitimate use of "stop" is "then stop", but it has to be in the non-initial schedule


        """

        bad_events = []
        for event_num, event in enumerate(schedule.events):
            if event_num > 0:
                event_strucs = event.get_recursive_componenets()
                if len(event_strucs) == 1:
                    struc = event_strucs[0]
                    label = struc.label
                    if label in ('SITE', 'ROUTE'):
                        prop = label.lower()
                        if not schedule.events[0].get_property_values(prop):
                            schedule.events[0].add_property_value(prop, struc)
                        else:
                            event.remove_property(prop, struc)
                        bad_events.append(event)
                if sched_num == 0 and event.directive and event.directive.value == 'stop' and event not in bad_events:
                    bad_events.append(event)
        for event in bad_events:
            schedule.remove_event(event)

    def modify_maxdose(schedule):
        """ Insert maxdose.directive if missing, for better matches with dictionary.
        """

        if schedule.maxdose and not schedule.maxdose.directive:
            if schedule.instruction.primary_directive:
                schedule.maxdose.directive = schedule.instruction.primary_directive
            elif schedule.maxdose.form:
                directive_name = Form.value_to_likely_directive.get(form.value, '')
                if directive_name:
                    schedule.maxdose.directive = Directive(directive_name)

    def modify_taper(schedule):
        """ Insert form into taper if missing. E.g.
            "ne at bedtime increase by 1 every day to reach effect or max 10 capsules" where we only can deduce that we increase by 1 capsule from the max expression
        """

        if schedule.taper and not schedule.taper.dose.form and schedule.instruction.primary_form:
            form_value = schedule.instruction.primary_form.value
            if schedule.taper.dose.quant.value > 1:
                plurality = 'plural'
            else:
                plurality = 'singular'
            form = Form(form_value, plurality, [])
            schedule.taper.dose.form = form


    def modify_dose_plurality(event):
        if event.dose:
            quant = event.dose.quant
            form = event.dose.form
            if quant.num_type in ('int', 'frac', 'decimal'):
                if quant.value <= 1:
                    form.plurality = 'singular'
                else:
                    form.plurality = 'plural'
            elif quant.num_type == 'range':
                if quant.low.num_type == 'var':
                    # For dictionary atoms, the plurality of the form is whatever it is in the string.
                    return
                form.plurality = 'plural'
            elif quant.num_type == 'var':
                # For dictionary atoms, the plurality of the form is whatever it is in the string.
                return

    def remove_repetitive_freq_and_periodicity(schedule):
        """ We sometimes have repetitive "daily" in the same schedule. E.g. "take once daily every day". Remove identical ones.
        """

        def remove_duplicates(alist):
            if not len(alist) > 1:
                return
            struc2representation = dict([(struc, struc.pprint()) for struc in alist])
            if len(struc2representation.values()) > len(set(struc2representation.values())):
                unique_reps = []
                remove_strucs = []
                for i, struc in enumerate(alist):
                    representation = struc2representation[struc]
                    if representation in unique_reps:
                        remove_strucs.append(struc)
                    else:
                        unique_reps.append(representation)

                for struc in remove_strucs:
                    alist.remove(struc)
                    struc.accounted_for_by_sem = None
                    struc.rules_used.append('*removed_from_sem_in_remove_repetitive_freq_and_periodicity*')

        remove_duplicates(schedule.periodicity)
        remove_duplicates(schedule.freq)
        return

    def modify_multiple_durations(sched_num, schedule):
        """ Modify the second Schedule with Duration = N days to be "For the next N Days" if in the previous schedule there was a calendar_event.
            Motivation:
            In "Take 2 tabs today, then take 1 tab for 4 days" the second schedule really means "take 1 tab for the NEXT 4 days."
        """

        if sched_num > 0 and schedule.duration and not schedule.duration.offset and schedule.duration.on_off is None and not schedule.duration.up_to:
            prev_schedule = schedule.instruction.schedules[sched_num - 1]
            time_unit_value = schedule.duration.time_interval.time_unit.value
            if not prev_schedule.calendar_event or prev_schedule.calendar_event[0].typ not in ('numeric', 'relative'):
                return
            if prev_schedule.calendar_event[0].typ == 'numeric' and prev_schedule.calendar_event[0].time_unit.value != time_unit_value:
                return
            schedule.duration.offset = 'next'
        return

    def modify_primary_directive_for_on_off_duration(schedule):
        """ Special case: in "on for 12 hours / off for 12 hours" case make sure directive = 'apply' and form = 'patch' in order to match the atoms
        """

        if (schedule.duration and schedule.duration.on_off and
            (not schedule.instruction.primary_directive or schedule.instruction.primary_directive.value in ('use', 'remove')) and
            (not schedule.instruction.primary_form or schedule.instruction.primary_form.value == 'patch')):
            schedule.instruction.primary_directive = Directive('apply')
            schedule.instruction.primary_directive.rules_used.append('*deduced_in_modify_primary_directive_for_on_off_duration*')
            if not schedule.instruction.primary_form:
               schedule.instruction.primary_form = Form(form_name = 'patch', plurality = '', constituents = [])
               schedule.instruction.primary_form.rules_used.append('*deduced_in_modify_primary_directive_for_on_off_duration*')
        return


    for parse in sig.parses:
        for instruction in parse.instructions:
            if not isinstance(instruction, DrugAdmin):
                continue
            determine_primary_directive_and_form(instruction)
            remove_duplicate_indications(instruction)
            #remove_extra_schedules(instruction)
            for sched_num, schedule in enumerate(instruction.schedules):
                modify_maxdose(schedule)
                modify_taper(schedule)
                modify_multiple_durations(sched_num, schedule)
                remove_repetitive_freq_and_periodicity(schedule)
                modify_primary_directive_for_on_off_duration(schedule)
                remove_specious_events(sched_num, schedule)
                for event in schedule.events:
                    modify_directive_from_route(event.directive, event.route, instruction.primary_form)
                    modify_dose_plurality(event)
                    insert_directive_into_event(event)
                    insert_substrate_into_event(event)
                    modify_nightly_bedtime_timings(event)
                    modify_specific_day(event)
                    modify_directive_in_event(event)

        parse_flattened = parse.flatten()
        if 'QUANT' in parse_flattened:
            insert_dose_into_event_from_quant(parse)
        if 'AND_CONJ' in parse_flattened:
            remove_and_conj_between_indications(parse)
        if 'max' in parse_flattened:
            remove_strucs_following_unparsed_maxdose(parse)

    return

################## END modify_sems()

################## START MATCHING PARSE TO DICTIONARY ATOMS

def reorder_strucs_in_canonical_order_and_omit_punctuation(struc_list):
    """ Reorders a list of Strucs in the canonical order specified in the global constant canonical_struc_order. Removes all spaces, and returns a new list.

        To take care of multiple strucs with the same label (e.g. Timing, Freq, Periodicity) we use additional sorting parameters: get_key(struc) and string representation.
    """

    labeled_struc_list = [struc for struc in struc_list if not struc.is_space_or_punctuation_only()]
    labeled_struc_list.sort(key = lambda struc: (canonical_struc_order.index(struc.label) if struc.label in canonical_struc_order else 100, get_key(struc), struc.quick_print_struc()))
    return labeled_struc_list

def get_key(struc):
    """ Creates a labeled string key for a top-level struc. E.g., for FREQ would generate "FREQ/range-day" if FREQ was of the form "2-4 times a day".

    It adds the label only at the top level.
    """

    if not struc:
        return None

    key = struc.string_key()
    if not key:
        return None
    if key is True:
        key = struc.label
    else:
        key = struc.label + '/' + key
    return key

def get_key_set(struc_list):
    """ Uses struc.string_key(self) to generate a set of string keys (one key per struc) that represent the maximal information that is minimally
        necessary for a sig key to match to a dictionary  entry.

        The key_set of each atom compatible (i.e that could be used for transducing) a given sig should be a subset of the key_set of the sig.
    """

    keys = set()
    for struc in struc_list:
        key = get_key(struc)
        if key:
            keys.add(key)

    return keys

def get_normalized_key_set(struc_list):
    """ For debugging/presentation only. Sorts the key set by the canonical order of strucs. Returns a string representation.
    """


    sem_strucs = [struc for struc in struc_list if struc.accounted_for_by_sem]
    struc_list = reorder_strucs_in_canonical_order_and_omit_punctuation(sem_strucs)

    keys = []
    for struc in struc_list:
        key = get_key(struc)
        if key and key not in keys:
            keys.append(key)

    return ' '.join(keys)

def get_compatible_atoms(struc_list):
    """ Takes a list of Strucs and returns the list of all the atoms that are semantically compatible with the list, to minimize the search space.

        Motivation for key matching:
        -   If a dict atom has, for example, a FREQ component and no part of the sig's SEM has FREQ, that atom can never be part of any version transducing the sig.
        -   If atom has FREQ/range-week and no FREQ segment of the sig's SEM has both range and week in FREQ, then that atom is out of consideration.
        -   If atom has directive 'use', it can stay in cosideration for being used for the version. But if the atom has directive 'instill', it could potentially be used only
            if the sig has a primary_directive of 'instill' OR 'use' (e.g. "use 1 vial" in a sig could match the atom "instill 1 via via nebulizer").

    We use the attribute sig.key_set that is attached to all atoms in the dictionary.
    """

    key_set = get_key_set(struc_list)
    all_atoms = atoms_dictionary.sigs
    compatible_atoms = [atom for atom in all_atoms if atom.key_set <= key_set]
    compatible_atoms = [atom for atom in compatible_atoms if atom.parses[0].is_fully_parsed()]

    if debug:
        #### Should not be needed by production time because we will need to have finished parsing all dictionary atoms
        compatible_atoms = [atom for atom in compatible_atoms if atom.parses[0].is_fully_parsed()]


    if debug:
        struc_rep = [struc.pprint() for struc in struc_list]
        error_msg = ('Found %d atoms in dictionary with compatible key_set to -->%s<--' % (len(compatible_atoms), struc_rep))
        msg = DeveloperMessages(error_msg, msg_type = 'Info', error_level = 0)
        errors_list.append(msg)


    return compatible_atoms


def prune_out_incompatible_atoms_by_keys(a_parse):
    """ Given a parse of a sig, prune out of the atomic dictionary all the atoms that can't possibly be of relevance, to minimize the search space. Returns a list of the compatible atoms.

        Motivation for key matching:
        -   If a dict atom has, for example, a FREQ component and no part of the sig's SEM has FREQ, that atom can never be part of any version transducing the sig.
        -   If atom has FREQ/range-week and no FREQ segment of the sig's SEM has both range and week in FREQ, then that atom is out of consideration.
        -   If atom has directive 'use', it can stay in cosideration for being used for the version. But if the atom has directive 'instill', it could potentially be used only
            if the sig has a primary_directive of 'instill' OR 'use' (e.g. "use 1 vial" in a sig could match the atom "instill 1 via via nebulizer").
    """

    struc_list = a_parse.strucs[:]
    labels = get_struc_labels(struc_list, delimiter = None, omit_spaces_and_punctuation = True)
    if 'DIRECTIVE' not in labels:
        for instruction in a_parse.instructions:
            directive = instruction.primary_directive
            if directive:
                struc_list.append(directive)
                break
    potentially_compatible_atoms = get_compatible_atoms(struc_list)
    return potentially_compatible_atoms

def get_sorted_compatible_atoms(a_parse):
    """ For debugging/presentation only.
    """

    atoms = prune_out_incompatible_atoms_by_keys(a_parse)
    atoms = [atom for atom in atoms if atom.parses[0].is_fully_parsed()]
    atoms.sort(key = lambda x: x.raw_sig)
    return atoms


class InstAtom(object):
    """ Instantiated Atoms -- i.e. an atom from the dictionary, together with the strucs from the sig it instantiates and numerical maps that map
        variables such as <<NUM_2>> to values found in the corresponding struc, such as "2.5".

    atom                    Sig struc that is the atom from the Dictionary

    matching_strucs         List of strucs from the sig that have been matched (covered) by the strucs in this atom.

    quality_match_score     Integer representing the quality of the match of the strucs in the matching_strucs to the corresponding struc in the atom, with values ranging from 0 to 4.

    numeric_mappings        List of pairs  [(var1, string repr of value1), (var2, value2), ...]. E.g., ('<<NUM_1>>', '1.5'), or ('<<DATE_0>>', '2/3/2012'),
                            where each var is taken from a struc in atom with a Quant sub-property.

    rank                    An integer that represents how important this atom is compared to other atoms in a version. It is higher if the atom has a DOSE struc
                            and if it is longer (has more strucs).

    core_matched_strucs     Set of sig strucs from matching_strucs that can't be repeated in mulitple atoms for the same version (i.e. self.matching_strucs minus Directive and Anaphora)


    """
    def __init__(self, atom, matching_strucs, quality_match_score, numeric_mappings):
        self.atom = atom
        self.matching_strucs = matching_strucs
        self.quality_match_score = quality_match_score
        self.numeric_mappings = numeric_mappings
        self.rank = 0
        self.core_matched_strucs = self.get_core_matched_strucs()

    def get_core_matched_strucs(self):
        return set([struc for struc in self.matching_strucs if struc.label not in ('DIRECTIVE', 'ANAPHORA')])

    def get_raw_rank(self):
        """ In the order of Atoms in the final version of transduction of an Event, the following are priorities, in descending order:
            -   Any atom that starts with AND_CONJ ("Also") or TAPER should come first
            -   Anything with DOSE should come first.
            -   The longer ones should come first.
            -   The ones first in canonical_struc_order should come before those later.
            -   Two Timing ones should be in the order they were encountered in the original text.
        """

        and_conj_rank = 0
        dose_rank = 0
        timing_pos_rank = 0
        position_rank = 0
        len_rank = len(self.core_matched_strucs)
        for struc in self.matching_strucs:
            label = struc.label
            if label == 'DOSE':
                dose_rank = 1
            elif label in ('AND_CONJ', 'TAPER'):
                and_conj_rank = 1
            elif label == 'AS_DIRECTED' and len_rank > 1:
                # AS_Directed combined with another Struc should have the same rank as that struc.
                # E.g. "USE AS DIRECTED EVERY DAY." should go after "USE THIS MEDICINE 2 TIMES A DAY."
                len_rank -= 1
            elif label == 'TIMING' and struc.accounted_for_by_sem:
                timing_pos_rank = struc.accounted_for_by_sem.timing.index(struc)

            if label not in ('ANAPHORA', 'DIRECTIVE'):
                position_rank += canonical_struc_order.index(label)


        score = 20000*and_conj_rank + 10000*dose_rank + 100*len_rank - position_rank - timing_pos_rank
        return score

    def is_incompatible_with(self, other_instatom):
        """ Return True if the two instatoms are incompatible, i.e. they share a non-core structure.

        E.g. DIRECTIVE DOSE ROUTE FREQ  is incompatible with DIRECTIVE DOSE PERIODICITY but is compatible with DIRETIVE ANAPHORA PERIODICITY
        """

        if self.core_matched_strucs & other_instatom.core_matched_strucs:
            return True
        else:
            return False

    def get_compatible_instatoms(self, instatoms_list):
        """ compatibles are instatoms that don't haven any core_matched_strucs in common, i.e. the pairwise intersection of self.core_matched_strucs is empty.
        """
        return [instatom for instatom in instatoms_list if not self.is_incompatible_with(instatom)]

    def instantiate_atom_for_locale(self, locale):
        """ Given a locale, returns the translation string of the atom for the locale, with numerical variables replaced by their instantiated values from the sig.
        """

        if locale not in self.atom.dictionary:
            return ''
        translation_with_variables = self.atom.dictionary[locale]
        atom_translation = instantiate_atom_string_numerically(translation_with_variables, self.numeric_mappings, locale)
        return atom_translation

    def get_raw_atom_string(self):
        """ Used in debugging and quality comparisons         """

        return self.atom.raw_sig

    def show(self):
        atom_descr = ('\nAtom %s' % self.get_raw_atom_string())
        atom_descr += ('\n     Quality: %s' % (self.quality_match_score))
        atom_descr += ('\n     Key Set: %s' % get_normalized_key_set(self.atom.parses[0].strucs))
        atom_descr += ('\n     Rank: %s' % self.get_raw_rank())
        numeric_mappings =  self.numeric_mappings
        atom_instantiated = self.instantiate_atom_for_locale(en)
        atom_descr += ('\n     Instantiated Atom: %s' % atom_instantiated)
        atom_instantiated_russian = self.instantiate_atom_for_locale(ru)
        atom_descr += ('\n     Instantiated Russian: %s' % atom_instantiated_russian)
        return atom_descr

    def pprint(self):
        print self.show()

class Version(object):
    """ Version of transduction of the Sig into a collection of dictionary Atoms.
        Represents one possible result of transduction.

    instatoms               List of InstAtoms (instantiated Atoms from Dictionary) that make up this transduction version

    quality_penalty         A non-negative integer that can be used to modify the quality score of the version.

    quality_scores_list     List of integers (for debugging): (quality_of_atoms_penalty, unmatched_strucs_penality, number_of_atoms_penalty, self.quality_penalty)
                            It is calculated by the self.quality_score() function.

    matched_strucs          Set of strucs from the sig's parse that have been matched  (covered) by at least one atom in this version.
                            May include some strucs that are not directly in parse.strucs list but that were deduced to be implied by it's semantics.

    """

    def __init__(self, parse):
        self.parse = parse
        self.instatoms = []
        self.quality_penalty = 0
        self.quality_scores_list = ()
        self.matched_strucs = set()

    def __add__(self, other):
        """ Implements concatentation of 2 versions: ver_3 = ver1 + ver2.  Returns a new object.
        """
        result = Version(self.parse)
        result.instatoms = self.instatoms + other.instatoms
        result.quality_penalty = self.quality_penalty + other.quality_penalty
        result.matched_strucs = self.matched_strucs | other.matched_strucs
        return result

    def copy(self):
        new_version = Version(self.parse)
        new_version.instatoms = self.instatoms[:]
        new_version.quality_penalty = self.quality_penalty
        new_version.matched_strucs |= self.matched_strucs
        return new_version

    def is_empty(self):
        return False if self.instatoms else True


    def add_instatoms(self, instatoms):
        """ Creates a new version from this version by adding each of the InstAtom (Instantiated Atom) objects
            from the instatoms and returning the new version. Except when self is a blank version, then it adds the instatoms and returns self.
        """

        new_version= self.copy()
        new_version.instatoms.extend(instatoms)
        new_version.matched_strucs.update(*[set(instatom.matching_strucs) for instatom in instatoms])
        return new_version

    def reorder_instatoms_in_schedule(self):
        """ Assuming this Version has only instatoms from one schedule, reorders them.

        Justification: need to pull upfront atoms such as:
        -   THEN_CHRONO CALENDAR_EVENT: (e.g., THEN ON DAYS <<NUM_0>>-<<NUM_1>>:)
        -   CALENDAR_EVENT: (e.g. ON DAYS <<NUM_0>>-<<NUM_1>>:)

        TODO: deal with atom "IF NO RELIEF:" when struc is created.
        """

        move_upfront = None
        for i, instatom in enumerate(self.instatoms):
            parse_string = instatom.atom.parses[0].flatten()
            if 'CALENDAR_EVENT:' in parse_string or 'DURATION:' in parse_string or 'THEN_CHRONO' in parse_string:
                move_upfront = instatom
                break
        if move_upfront and i > 0:
            del self.instatoms[i:i+1]
            self.instatoms.insert(0, move_upfront)

    def reorder_instatoms_in_instruction(self):
        """ Assuming this Version has only instatoms from one instruction, reorders them.

        Justification:
        -   REPEAT and MISC (type = call_911) mention pain, so the Indication needs to move before them.

        """

        indication_pos = None
        repeat_or_misc_pos = None
        for i, instatom in enumerate(self.instatoms):
            parse_string = instatom.atom.parses[0].flatten()
            if 'REPEAT' in parse_string or 'MISCELLANEOUS' in parse_string and repeat_or_misc_pos is None:
                repeat_or_misc_pos = i
            if 'INDICATION' in parse_string and indication_pos is None:
                indication_pos = i
                indication = instatom

        if indication_pos is not None and  repeat_or_misc_pos is not None and repeat_or_misc_pos < indication_pos:
            del self.instatoms[indication_pos:indication_pos+1]
            self.instatoms.insert(repeat_or_misc_pos, indication)


    def cleanup_remove_debris(self):
        """ Remove terminals "AFTER THAT:" and "ALSO:"
        """

        remove_after = None
        num_instatoms = len(self.instatoms)
        for i in range(num_instatoms):
            last_phrase = self.instatoms[-1 - i].atom.raw_sig
            if last_phrase in ('AFTER THAT:', 'ALSO:'):
                remove_after = num_instatoms - 1 - i
            else:
                break
        if remove_after is not None:
            self.instatoms[remove_after:] = []

    def get_unmatched_strucs(self, struc_list):
        unmatched_strucs = [struc for struc in struc_list if struc not in self.matched_strucs]
        return unmatched_strucs

    def strucs_not_covered(self):
        """ Returns the list of strucs in the Sem that are not covered by this version.
        """

        all_strucs = self.parse.get_sem_componenets()
        unmatched_strucs = self.get_unmatched_strucs(all_strucs)

        return unmatched_strucs

    def assess_quality_penalty_for_event(self, event, is_initial_event):
        """ Called when the version represents (at this point of processing) only one event (though it can have Schedule strucs, too, if there is 1 event in the Schedule).

            This is needed because sometimes you can't rearrange instatoms to make the version work.
            E.g. if instatom has both "Also" and TIMING but another instatom has DOSE, no arrangement
            will make it right.

            Ensures, for example, that a version that has Dose instatom after Timing is penalized.

            One exception to this is the non-initial doses. There, timing should come first (ideally with the Dose, but not after it).
            E.g. "take 3 capsules in the morning and 3 capsules at 6pm" should NOT transduce to "take 3 capsules in the morning. Also: take 3 capsules. Take this medicien at 18:00",
            Rather, it shloud transduce to "take 3 capsules in the morning. Also: take this medicine at 18:00. Take 3 capsules"

            Other reward/penalties:
            -   The longer the instatom with DOSE, the better.

            -   Substituting "affected site" for actual site should lower the score.

            -   It is good to have "SITE" with DOSE.  Sounds better.

            -   Penalize those version where the first atom does not accurately reflect the directive of the sig (e.g. uses the more general "use" vs. specific "inject").
                Although exact match is always rewarded automatically in directive.match_dictionary(), it is much more important that the match be exact in the first atom.
                Example: "swish and spit four times daily": if we have a directive "swish and spit" anywhere in the dictionary, it is important to use it in the first atom.
                But is ok if "use every day" is used in transduction instead of "swish and spit every day"

            -   DOSE should not be in the same sentence as PERIODICITY = day if there is FREQ with time inteval = day.
                E.g. "take 2 tablets 3 times a day every day" can't have in its transduction an atom "take 2 tablets every day", "take this medicine 3 times a day"
                because in English "take 2 tablets every day" implies that 2 is the total daily dose, while in reality the total daily dose is 6.

            -   More generally, if there are several FREQ or PERIODICITY (e.g 2 FREQs or 1 FREQ and 1 PERIODICITY) in the same Schedule, then only the ones with the smallest time_unit can be together with the DOSE.
                - Example 1:  "take 2 tablets 3 times a day every day" can't have in its transduction an atom "take 2 tablets every day", "take this medicine 3 times a day"
                because in English "take 2 tablets every day" implies that 2 is the total daily dose, while in reality the total daily dose is 6.
                - Example 2:   "take 2 tabs daily every 4 hours") the lowest-period periodicity should be the one with the DOSE (if any is). I.e. the above sentence should NOT transduce to
                "take 2 tablets every day. Take this medicine every 4 hours' but to "take 2 tablets every 4 hours. take this medicine every day."
        """

        labels_covered_so_far = set()
        labels_that_should_not_preced_dose = set(['TIMING', 'VEHICLE', 'SITE', 'ROUTE'])
        dose_found = False

        for instatom_index, instatom in enumerate(self.instatoms):
            labels = [struc.label for struc in instatom.matching_strucs]
            if 'DOSE' in labels and not dose_found:
                if is_initial_event and not labels_covered_so_far.isdisjoint(labels_that_should_not_preced_dose):
                    self.quality_penalty += 3
                elif not is_initial_event and 'TIMING' in labels:
                    self.quality_penalty -= 2
                elif not is_initial_event and not 'TIMING' in labels_covered_so_far and event.timing:
                    self.quality_penalty += 3

                if 'PERIODICITY' in labels or 'FREQ' in labels:
                    # See "Other reward/penalties" above. We need to make sure that DOSE can only be in the phrase with the lowest time_unit Periodicity or FREQ.
                    # E.g. DOSE should not be in the same sentence as PERIODICITY = day if there is FREQ with time inteval = day.
                    instatom_pers_or_freqs = [struc for struc in instatom.matching_strucs if struc.label in ('PERIODICITY', 'FREQ')]
                    other_pers_or_freqs = [struc for struc in self.matched_strucs if struc.label in ('PERIODICITY', 'FREQ') and struc not in instatom_pers_or_freqs]
                    # make sure that a periodicity structure in this instatom has the period (time_unit) less than any other freq or periodicity not in the instatom.
                    # make sure that a freq structure in this instatom has the period (time_unit) less than or equal to any other freq or periodicity not in the instatom.
                    # E.g. "take 2 tabs 3 times daily" is ok if "every other day" is in another instatom. But "take 2 tabs every day" is NOT ok if "3 times a day" is in another instatom.

                    # TimeUnit.permissible_values are listed from smallest to largest: (second, minute, hour, day, week, month). So Dose can only be with the smallest such periods.

                    for local_struc in instatom_pers_or_freqs:
                        for other_struc in other_pers_or_freqs:
                            local_time_unit_index = TimeUnit.permissible_values.index(local_struc.time_unit.value)
                            other_struc_time_unit_index = TimeUnit.permissible_values.index(other_struc.time_unit.value)
                            if local_struc.label == 'PERIODICITY' and local_time_unit_index >= other_struc_time_unit_index:
                                if other_struc.quant.value == 1 and local_time_unit_index == other_struc_time_unit_index:
                                    # It is ok (though not great) to have Periodicity with the Dose if the other Freq (outside of instatom) has period = 1. E.g. "take 2 tablets every day. take this medicine once a day"
                                    # At least this transduction does not distort that actual daily dose one needs to take.
                                    self.quality_penalty += 1
                                else:
                                    self.quality_penalty += 6
                                    break
                            elif local_struc.label == 'FREQ' and local_time_unit_index > other_struc_time_unit_index:
                                self.quality_penalty += 6
                                break

                if 'SITE' in labels:
                    self.quality_penalty -= 0.5

                if instatom.quality_match_score >= 3.0:
                    # only reward the Dose instatom for length if it is of high quality
                    length_reward = len(labels) - 2
                else:
                    length_reward = 0
                self.quality_penalty -= length_reward
                dose_found = True

            if 'SITE' in labels:
                sig_site = [struc for struc in instatom.matching_strucs if struc.label == 'SITE'][0]
                instatom_site =  [struc for struc in instatom.atom.parses[0].strucs if struc.label == 'SITE']
                if not instatom_site or instatom_site[0].value != sig_site.value:
                    self.quality_penalty += 2

            if instatom_index == 0 and 'DIRECTIVE' in labels:
                # Penalize those version where the first atom does not accurately reflect the directive of the sig (e.g. "use" vs. "inject")
                # Although exact match is always rewarded automatically in directive.match_dictionary(), it is much more important that the match be exact in the first atom.
                sig_directive = [struc for struc in instatom.matching_strucs if struc.label == 'DIRECTIVE'][0]
                instatom_directives = [struc for struc in instatom.atom.parses[0].strucs if struc.label == 'DIRECTIVE']
                if instatom_directives:
                    instatom_directive = instatom_directives[0]
                    if sig_directive.value != instatom_directive.value:
                        self.quality_penalty += 0.5

            labels_covered_so_far |= set(labels)


    def quality_score(self):
        if self.is_empty():
            return 0
        else:
            #average_quality_of_atoms = round(sum([instatom.quality_match_score for instatom in self.instatoms])* 1.0 / len(self.instatoms), 2)
            #quality_of_atoms_penalty = 4 - average_quality_of_atoms
            quality_of_atoms_penalty = sum([4 - instatom.quality_match_score for instatom in self.instatoms]) * 0.5

            # the following penalty is 3 points for each struc in parse.strucs that has not been matched by any of the atoms
            # We except Then_Chrono and And_Cong because they are the only strucs in parse.strucs that are not explicitly
            # represented in the Sem representation of the Sig.

            unmatched_strucs_penality = 6 * len(self.strucs_not_covered())
            #unmatched_strucs_penality = 3 * len([struc for struc in self.parse.strucs if struc not in self.matched_strucs
            #                                     and struc.accounted_for_by_sem and struc.label not in ('THEN_CHRONO', 'AND_CONJ')])

            minimum_expected_number_of_atoms = sum([instruction.minimum_expected_number_of_atoms() for instruction in self.parse.instructions])
            number_of_atoms_penalty = 0.3 * (len(self.instatoms) - minimum_expected_number_of_atoms)
            self.quality_scores_list = (quality_of_atoms_penalty, unmatched_strucs_penality, number_of_atoms_penalty, self.quality_penalty)
            value = round(10 - quality_of_atoms_penalty - unmatched_strucs_penality - number_of_atoms_penalty - self.quality_penalty, 2)
            return value


    def get_raw_atoms_strings(self):
        """ Returns list of strings or raw atoms that comprize the version. Used for debugging and quality comparison purposes        """

        list_of_atom_strings = [instatom.get_raw_atom_string() for instatom in self.instatoms]
        return list_of_atom_strings

    def get_instantiated_atoms_list_for_locale(self, locale):
        """ Returns list of instantiated translation strings for a locale. E.g., for locale == en, may return ["Take 2-4 tablets 2 times a day.", "Take this medicine for acne."]

            Used with locale == en to find equivalent versions up to permutations of atoms.
        """

        list_of_instantiated_translation_strings = [instatom.instantiate_atom_for_locale(locale) for instatom in self.instatoms]
        return list_of_instantiated_translation_strings

    def get_instantiated_translations(self):
        """ Returns a dictionary mapping locales to translation string that represents the concatenated list of all atoms
        translated and numerically instantiated for the locale.        """

        locale_2_translation_string = dict()
        for locale in available_locales:
            list_of_instantiated_translation_strings = self.get_instantiated_atoms_list_for_locale(locale)
            translation = ' '.join(list_of_instantiated_translation_strings)
            locale_2_translation_string[locale] = translation

        return locale_2_translation_string

    def canonical_string_representation(self):
        """ Provides canonical string representaiton for Version, so that 2 Versions with the same canonical rep are substantially identical
            if they have the same string rep.

            The representation consists of the first atom's english version concatenated with the ordered list of the remaining atoms' Enlish version
        """
        strings_list = self.get_instantiated_atoms_list_for_locale(locale = en)
        if not strings_list:
            return ''
        header_atom = strings_list[0]
        tail_atoms = sorted(strings_list[1:])
        if tail_atoms:
            return header_atom + '|' + '|'.join(tail_atoms)
        else:
            return header_atom

    def __eq__(self, other):
        """ Two Versions are substantially the same if they instantiate to the same set of english atoms and the first atoms are identical.        """

        self_repr = self.canonical_string_representation()
        other_repr = other.canonical_string_representation()
        return self_repr == other_repr


    def __ne__(self, other):
        return not self.__eq__(other)



    @classmethod
    def concatenate_versions(cls, *version_lists):
        """ Creates a modified cross product of all version_lists by concatenating one version from each list pointwise.

            If some of the lists are empty, they are simply omitted from the product (a strict cross-product would have
            to return an empty list if one of the components was empty).
        """

        version_lists = [version_list for version_list in version_lists if version_list]

        if not version_lists:
            return version_lists
        elif len(version_lists) == 1:
            return version_lists[0]

        cross_product = cross_lists(*version_lists)

        parse = version_lists[0][0].parse
        new_versions = []
        for row in cross_product:
            new_version = Version(parse)
            for version in row:
                if not version.is_empty():
                    new_version = new_version + version
            new_versions.append(new_version)
        return new_versions



    @classmethod
    def prune_and_sort_versions(cls, versions_list):
        """ Prunes out substantially identical versions and returns a list of unique versions in the order of decreasing quality_score
        """

        versions_list.sort(key = lambda version: version.quality_score(), reverse = True)
        version_representations = []
        unique_versions = []
        for version in versions_list:
            version_repr = version.canonical_string_representation()
            if version_repr in version_representations:
                continue
            else:
                unique_versions.append(version)
                version_representations.append(version_repr)
        return unique_versions

    def show(self, ver_num = None):

        descr = ''
        num_of_atoms_in_version = len(self.instatoms)
        descr += ('--> Version%s: quality %s, number of atoms: %d' % ('' if ver_num is None else ' ' + str(ver_num),  self.quality_score(), num_of_atoms_in_version))
        all_atoms_instantiated = []
        for atom_num, instatom in enumerate(self.instatoms):
            atom_descr = ('\nAtom %d: %s%s' % (atom_num, instatom.atom.raw_sig, '' if len(instatom.atom.dictionary) > 1 else ' [-->From TENTATIVE DICT]'))
            atom_descr += ('\n     Quality: %s' % (instatom.quality_match_score))
            atom_descr += ('\n     Key Set: %s' % get_normalized_key_set(instatom.atom.parses[0].strucs))
            numeric_mappings =  instatom.numeric_mappings
            atom_instantiated = instatom.instantiate_atom_for_locale(en)
            atom_descr += ('\n     Instantiated Atom: %s' % atom_instantiated)
            all_atoms_instantiated.append(atom_instantiated)
            atom_instantiated_russian = instatom.instantiate_atom_for_locale(ru)
            atom_descr += ('\n     Instantiated Russian: %s' % atom_instantiated_russian)
            descr += padlines(atom_descr, 5)
        descr += padlines(('\n----\nVersion%s Final: %s' % ('' if ver_num is None else ' ' + str(ver_num), '|'.join(all_atoms_instantiated))), 5)
        descr += padlines(('\n%d strucs in Sig not covered by the version: %s' % (len(self.strucs_not_covered()), '|'.join([struc.quick_print_struc() for struc in self.strucs_not_covered()]))), 5)
        descr += padlines(('\nQuality penalty subscores: atoms: %s, unmatched strucs: %s, Num of atoms: %s, outside penalty: %s' % self.quality_scores_list), 5)
        matched_struc_list = reorder_strucs_in_canonical_order_and_omit_punctuation(self.matched_strucs)
        matched_strucs_descriptor = '\n'.join([struc.pprint() + '{' + struc.show_accounted_for_by_sem(omit_coords = False) + '}' for struc in matched_struc_list])
        descr += padlines(('\nMatched strucs: \n%s\n--------\n' % (padlines(matched_strucs_descriptor, 5))), 5)
        return descr

    def pprint(self, ver_num = None):
        print self.show(ver_num = ver_num)

    def __str__(self):
        self.pprint()




def break_sig_into_atoms(sig):
    """ Breaks up semmed sig into potential atoms. Creates in sig.versions a list of Versions of how the sig can be transduced to lists of atoms.

    Takes each Instruction in the Sem representation and finds all Versions (combinations of atoms from the dictionary) so that each Version maximally covers the meaning of each
    Struc in the sig's Instruction.

    Approach and motivation:
    Ideally we want to find all partitions of the Strucs in the sig into groups of Strucs and then try to match each group to the dictionary atoms that have the same strucs.
    But we can't ry all partitions of the Strucs even for the core (The first event strucs plus the 1st schedule strucs) because of the exponential growth of the number of partitions:
    even if there are only 7 strucs, there are 877 partitions of the set of strucs into subsets, while for 8 strucs there are 4140 partitions one would have to test against the
    dictionary.

    Moreover, some Strucs (Directive, Anaphora) can repeat in multiple InstAtoms in one Version (e.g. "Take 1 tablet once daily" (DIRECTIVE DOSE FREQ PERIODICITY)
    is best transduced by a version with 2 InstAtoms which produce overlapping partitions of the Strucs (overlapping in DIRECTIVE):
    "TAKE 1 TABLET ONCE A DAY." (DIRECTIVE DOSE FREQ) and  "TAKE THIS MEDICINE EVERY DAY." (DIRECTIVE ANAPHORA PERIODICITY).

    So to minimize the search space, we do this:
    -   Break up the Sig into Sem chunks, and do matching only with semantically meaningful subsets (chunks) of strucs in the Sig. E.g., we try the first Event and the strucs of the corresponding Schedule
        together, but for the second Event or for strucs proper to Instrucion we do separate matching to atoms.
    -   For each chunk of Sig, we use our key procedure match_sem_chunk_to_dictionary() to prune out the list of atoms to just the ones that can possibly be used in any Version that transduces this chunk of Sig.
    -   Now that we have a list of InstAtoms each of which adequately matches some part of the meaning of the sig chunk, we use a recursive procedure group_matching_instatoms_into_versions()
        to group these Atoms into maximal subsets (Versions) so that the Instatoms in each version are consistent with each other and that you can't add any instatom to the Version that covers any
        non-covered Struc and yet remain consistent with the other InstAtoms of the version.


    Definitions:
    -   Version: an object that, by the end of this procedure, holds a maximal list of Atoms that together translate a portion of the Sig. It is
        "maximal" in the sense that we can't add any more atoms to it that would cover any more Strucs from the Sig that are not already covered.
        One can have multiple Versions that partition the Sem strucs in different ways, or where some Versions cover some concepts but omit others.

    -   InstAtom: an Atom each of whose structures corresponds to a Struc from the Sig and whose variables (NUM_N) are instantiated with the values from the Sig.

    -   Frame: All atoms whose structure is represented by the same list of Struc labels. E.g. "TAKE 1 OR 2 TABLETS AT BEDTIME." and "TAKE 1 CAPSULE IN THE MORNING." have the same
        frame because their strucs have the same sequence of labels: "DIRECTIVE DOSE TIMING".
    """

    debug_recursive_calls = [0]    # A facility to count the number of recursive calls (not exactly depth in our case) to find_maximal_instatom_groups_recursively()


    def get_structurally_similar_list_of_atoms_in_dict(struc_list):
        """ Get the list of atoms that have the same strucs in the same order as the given struc_list after being put in canonical order.
        """

        canonically_ordered_struc_list = reorder_strucs_in_canonical_order_and_omit_punctuation(struc_list)
        struc_list_string = get_struc_labels(canonically_ordered_struc_list, delimiter = '|', omit_spaces_and_punctuation = False)
        candidates = atoms_dictionary.match_labels_list_2_atoms_list.get(struc_list_string, [])
        if debug:
            error_msg = ('Found %d candidate atoms in dictionary matching struc description %s' % (len(candidates), struc_list_string))
            msg = DeveloperMessages(error_msg, msg_type = 'Info', error_level = 0)
            errors_list.append(msg)
        return candidates

    def match_struc_list_to_atom(struc_list, atom):
        """ Matches struc_list to the atom. Returns (match_quality, numeric_mappings) pair if match is good enough, otherwise returns None

            match_quality       Float between 4 and 0
            numeric_mappings    List of pairs (var name, string repr of value)
        """

        def debugging_find_substruc_with_lowest_quality_score(sig_struc, dict_struc):
            final_match_quality = sig_struc.match_dictionary(dict_struc)
            for prop in sig_struc.substantive_properties:
                if prop not in dict_struc.substantive_properties:
                    continue
                try:
                    sig_prop_value = sig_struc.__dict__[prop]
                    dict_prop_value = dict_struc.__dict__[prop]
                    match_quality = sig_prop_value.match_dictionary(dict_prop_value)
                    if match_quality == final_match_quality:
                        return (sig_prop_value, dict_prop_value)
                except AttributeError:
                    continue
            return None

        if not struc_list:
            return None
        atom_struc_list = atom.parses[0].strucs
        atom_strucs = reorder_strucs_in_canonical_order_and_omit_punctuation(atom_struc_list)
        labels_list = [struc.label for struc in struc_list]
        combined_match_quality = []
        numeric_mappings = []
        for struc_num, sig_struc in enumerate(struc_list):
            if struc_num >= len(atom_strucs):
                if debug:
                    raise Exception('Unexpected mismatch of atom length and struc list length in match_struc_list_to_atom()')
                return None
            atom_struc = atom_strucs[struc_num]
            if sig_struc.label != atom_struc.label:
                if debug:
                    raise Exception('Unexpected mismatch of atom labels and struc list labels in match_struc_list_to_atom()')
                return None

            match_quality = sig_struc.match_dictionary(atom_struc)

            if sig_struc.label == 'ANAPHORA' and len(struc_list) > 1:
                # Don't dilute the scores of other strucs by anaphora matches. Anaphora either matches or it doesn't.
                pass
            elif sig_struc.label == 'DIRECTIVE' and match_quality < 4 and match_quality > 1 and 'DOSE' not in labels_list:
                # This is probably "Use" instead of a more precise directive, or "take" istead of "give", which get quality score of 2.
                # It is more important that the precise directive be included with the DOSE instatom than with supplementary instatoms, so reduce
                # penalty in the latter case by half.
                modified_match_quality = 4 - 0.5 * (4 - match_quality)
                combined_match_quality.append(modified_match_quality)
            else:
                combined_match_quality.append(match_quality)

            if match_quality < 1:
                if debug:
                    error_msg = ('--> In match_struc_list_to_atom() with struc_list %s with atom %s. Bad quality for struc %s' %
                                 (get_struc_labels(struc_list, delimiter = ' ', omit_spaces_and_punctuation = True), atom.raw_sig, atom_struc.label))
                    result = debugging_find_substruc_with_lowest_quality_score(sig_struc, atom_struc)
                    if result:
                        (sig_prop_value, dict_prop_value) = result
                        error_msg += ('\n       Error is in substruc %s of %s.\n      Sig_Struc: \n%s\n     Atom_Struc:\n%s' %
                                      (sig_prop_value.label, atom_struc.label, padlines(sig_prop_value.__str__(), 10), padlines(dict_prop_value.__str__(), 10)))
                    msg = DeveloperMessages(error_msg, msg_type = 'Info', error_level = 0)
                    errors_list.append(msg)
                break
            numeric_mapping = sig_struc.get_numerical_map_to_dict(atom_struc)
            numeric_mappings += numeric_mapping

        if combined_match_quality and min(combined_match_quality) >= 1:
            #min_match_quality = min(combined_match_quality)
            #average_match_quality = sum(combined_match_quality) / float(len(combined_match_quality))
            #match_quality = round((min_match_quality + average_match_quality)/2, 2)
            match_quality_deficit = sum([4 - quality for quality in combined_match_quality])
            match_quality = 4 - match_quality_deficit
            return (match_quality, numeric_mappings)
        else:
            return None


    def try_to_match_tentative_struc_list_to_dict(tentative_struc_list, restrict_to_these_atoms_list = None, max_versions = 2):
        """ Tries to match a list of strucs to dictionary. Returns the list of lists: (atom, quality_score, numeric_mappings) with quality_score > 0
            sorted in decreasing order of quality.

        restrict_to_these_atoms_list    List of atoms to which the output is to be restricted

        max_versions                    INT: limits the number of atoms with identical struc that are returned.

        numeric_mappings                List of pairs from variable name in dictionary to it's numeric value as string.
        """

        if not restrict_to_these_atoms_list:
            restrict_to_these_atoms_list = potentially_applicable_atoms_for_whole_sig

        candidate_atoms_by_structure = get_structurally_similar_list_of_atoms_in_dict(tentative_struc_list)
        if len(candidate_atoms_by_structure) < len(restrict_to_these_atoms_list):
            candidate_atoms = [atom for atom in candidate_atoms_by_structure if atom in restrict_to_these_atoms_list]
        else:
            candidate_atoms = [atom for atom in restrict_to_these_atoms_list if atom in candidate_atoms_by_structure]

        if debug:
            error_msg = ('In try_to_match_tentative_struc_list_to_dict() Found %d candidate_atoms for struc_list %s' %
                         (len(candidate_atoms), get_struc_labels(tentative_struc_list, delimiter = ' ', omit_spaces_and_punctuation = True)))
            msg = DeveloperMessages(error_msg, msg_type = 'Info', error_level = 0)
            errors_list.append(msg)


        reordered_tentative_struc_list = reorder_strucs_in_canonical_order_and_omit_punctuation(tentative_struc_list)
        good_enough_matches = []

        for atom in candidate_atoms:
            match_result = match_struc_list_to_atom(reordered_tentative_struc_list, atom)
            if match_result:
                (match_quality, numeric_mappings) = match_result
                good_enough_matches.append((atom, match_quality, numeric_mappings))

        good_enough_matches.sort(key = lambda x: x[1], reverse = True)

        if debug:
            error_msg = ('--> Found %d atoms in dictionary with good enough quality match to %s' % (len(good_enough_matches), get_struc_labels(tentative_struc_list, delimiter = ' ', omit_spaces_and_punctuation = True)))
            msg = DeveloperMessages(error_msg, msg_type = 'Info', error_level = 0)
            errors_list.append(msg)
            #print error_msg
            for (atom, quality_score, numeric_mappings) in good_enough_matches:
                pass
                #print quality_score, atom.raw_sig, numeric_mappings

        good_enough_matches = good_enough_matches[:max_versions]

        return good_enough_matches

    def try_adding_struc_list_to_a_version(struc_list, version, restrict_to_these_atoms_list):
        """ Takes a list of strucs (some of which may be part of parse.strucs and some are new strucs) and tries to match it to the atoms dictionary.
            If successful, returns a list of new versions, one for each possible dictionary match.

            If some strucs in struc_list are None, we just omit them.

            If struc_list doesn't match any atoms well enough, return []
        """

        struc_list = [struc for struc in struc_list if struc]
        if not struc_list:
            return []

        candidate_atoms = try_to_match_tentative_struc_list_to_dict(struc_list, restrict_to_these_atoms_list)
        instatoms = [InstAtom(atom, struc_list, quality_match_score, numeric_mappings) for (atom, quality_match_score, numeric_mappings) in candidate_atoms]
        new_versions = [version.add_instatoms([instatom]) for instatom in instatoms]
        return new_versions


    def try_adding_struc_list_to_versions_list(struc_list, versions_list, restrict_to_these_atoms_list = None):
        """ Try to extend each version in the versions_list by trying to add to it an atom with Frame = struc_list. Return updated_versions.

            Return the list of new versions that are actual new extensions of some version from the version_list.
            If a version can't be extended, it is not part of the returned list.
        """

        extended_versions = []
        for version in versions_list:
            matched_strucs = version.matched_strucs
            if set(struc_list) <= matched_strucs:
                continue
            else:
                new_extended_versions = try_adding_struc_list_to_a_version(struc_list, version, restrict_to_these_atoms_list)
                extended_versions += new_extended_versions
        return extended_versions


    def find_portion_of_struc_list_corresponding_to_atom(struc_list, atom_strucs, atom_strucs_labels):
        """ Finds the portion of the struc list that corresponds label-for-label to the atom strucs. Returns the list of found matching_strucs (each
            a sublist of struc_list). (Usually returns just singleton list [matching_strucs], but cardinality could be greater if multiple matches are found).

            If a struc.label occurs multiple times in either atom or struc_list, the multiplicities have to match on both sides.
            If a struc.label occurs mulitple times in the struc_list but fewer times in the atom_strucs, we need to return all possible matches.
            Example 1: if atom is "APPLY EVERY OTHER NIGHT AT BEDTIME." with atom_strucs = DIRECTIVE PERIODICITY TIMING TIMING, then struc_list better have
            two Timing strucs.
            Example 2: But if atom is "TAKE 1 TABLET AT BEDTIME." (DIRECTIVE DOSE TIMING with only one Timing), while the struc_list is DIRECTIVE DOSE PERIODICITY TIMING TIMING
            ("take 1 tablet every night at bedtime"), then we need to output 2 possible matches DIRECTIVE DOSE TIMING for the atom_strucs:
            one where Timing = "at night", the other where Timing = 'bedtime'.
            Example 3: Multiple CALENDAR_EVENTs in the sig: e.g. "Take 2 tablets at once on Day 1 ...".
            Example 4: Three Timing strucs in the sig: "take one tablet by mouth in the morning on an empty stomach 30 minutes before food"

            Both struc_list and atom_strucs have been sorted via reorder_strucs_in_canonical_order_and_omit_punctuation() prior to being
            sent to this procedure.
        """

        # We need to assume that each struc.label appears only once in the struc_list. Otherwise we need to deal
        # with needing to reorder the multiple Timing or Freq structures to make them match the atom. I.e. we will need to deal with multiple Timings, Freqs, Periodicity, and
        # Calendar_Events in one event/schedule separately.
        matching_strucs = [struc for struc in struc_list if struc.label in atom_strucs_labels]
        matching_strucs_labels = [struc.label for struc in matching_strucs]
        if matching_strucs_labels ==  atom_strucs_labels:
            return [matching_strucs]
        if set(atom_strucs_labels) != set(matching_strucs_labels):  # This can happen is the atom has a struc not in sig, e.g. "and_conj"
            return []

        # We have a case of labels with multiplicity > 1 on one or both sides, where the multiplicities from the two sides don't match.
        # For each label, let it's multiplicity in struc_list = S and it's mulitplicity in atom = A:
        #   -   if S < A, no match is possible.
        #   -   if S > A, find all sublists of strucs with that label of length A.
        #   -   if S == A, the match is trivial.
        # So each actual matching_strucs list will be composed of concatenation of these sublists: a run of strucs for which labels have S==A, then one sublist from
        # the list of sublists of struc_list of length A for the label for which S > A, etc. I.e. we do cross_product from lists of sublists

        sig_label_2_multiplicity = defaultdict(int)
        for label in matching_strucs_labels:
            sig_label_2_multiplicity[label] += 1

        atom_label_2_multiplicity = defaultdict(int)
        for label in atom_strucs_labels:
            atom_label_2_multiplicity[label] += 1


        segments = []               # Each segment is a list. It is either list of length 1 and contains a list of all conseq strucs for which S == A, (segment = [[struc1, struc2, ...]])
                                    # Or it is a list of lists of length A (if S > A) of strucs with the same label (segment = [[struc1, struc2], [struc2, struc3], [struc1, struc3]] of A == 2 and S == 3)
        start_new_segment = True    # A flag that means the preceding struc has S > A and, if the current struc has S==A, we now need to start a new segment.
                                    # We reset the flag to True whenever we are on a structure with S > A.
        repeating_label = None      # Label for the case S > A
        for struc in matching_strucs:
            if sig_label_2_multiplicity[struc.label] < atom_label_2_multiplicity[struc.label]:
                return []
            if sig_label_2_multiplicity[struc.label] == atom_label_2_multiplicity[struc.label]:
                if repeating_label:     # We are starting a new segment with A == S, but the previous struc had S > A, so close off that segment
                    previous_segment = choose_m_of_n(segment, atom_label_2_multiplicity[repeating_label])
                    segments.append(previous_segment)
                    repeating_label = None
                if start_new_segment:
                    segment = [struc]
                    start_new_segment = False
                else:
                    segment.append(struc)
            else:                           # sig_label_2_multiplicity[struc.label] > atom_label_2_multiplicity[struc.label]
                if not start_new_segment:   # The previous struc had A==S and was in a segment that we now need to close and reset the flag to True
                    segments.append([segment])
                start_new_segment = True
                if struc.label == repeating_label:
                    segment.append(struc)
                    continue
                elif repeating_label:       # We are starting a new label for which S > A but the previous struc also had S > A, so close off that segment
                    previous_segment = choose_m_of_n(segment, atom_label_2_multiplicity[repeating_label])
                    segments.append(previous_segment)
                # We are starting a new label for which S > A and the previous struc had S == A
                repeating_label = struc.label
                segment = [struc]

        if not start_new_segment:   # The last struc in matching_strucs list had A == S, so close that last segment
            segments.append([segment])
        if repeating_label:         # The last struc had S > A, so close off that segment
            previous_segment = choose_m_of_n(segment, atom_label_2_multiplicity[repeating_label])
            segments.append(previous_segment)

        matching_strucs_list_of_lists = cross_lists(*segments)
        matching_strucs_list = []
        for alist in matching_strucs_list_of_lists:
            segment = []
            for struclist in alist:
                segment.extend(struclist)
            matching_strucs_list.append(segment)

        return matching_strucs_list



    def match_sem_chunk_to_dictionary(struc_list, instruction):
        """ Takes the struc_list and returns a list of Versions that transduce it.

        struc_list                  List of strucs from the sig. Typically may be all strucs related to an event that could be part of the transductions for the event,
                                    possibly minus some strucs that were previously removed, e.g. Duration or Calendar_Event if they were taken out to frame the event)
        instruction                 A reference to the instruction sem. Used only to get the reference to parse and to insert the right DIRECTIVE if needed.

        1.  Find all dictionary atoms compatible with struc_list (compatible_atoms).
        2.  From compatible_atoms find those atoms that match any portion of the struc_list or the anaphorized struc_list. When a transduction involves more than
            one atom, the subsequent atoms are likely to use anaphoras such as "Take THESE MEDICINES in the evening", etc. Moreover, you can't have
            an anaphora with a dose -- anaphora "this medication" typically refers to the dose from a previous phrase, so we remove DOSE struc from the anaphorized struc_list.
            Therefore, each atom that will be used in a transduction version should match exactly the frame consisting of strucs from the actual sig (struc_list) or the strucs
            from the anaphorized struc_list.
        3   Instantiate these compatible items to produce the list matched_atoms of Instatoms. These are all and only instatoms that can be used in various combinations to transduce the struc_list.
        4.  Prune matched_atoms so that if two atoms from struc_list match the same set of strucs, we keep only the one that matches best. We do this to keep the
            computational costs of the next step (the recursion to find maximal groups of these matched instatoms that are compatible with each other and
            provide maximal cover of strucs in the struc_list). The list of these instatoms is best_matched_atoms.
        5.  Group these matching Instatoms from best_matched_atoms into maximal groups maximal_instatom_groups. A maximal groups of these matched instatoms
            is a list of instatoms so that any two of them are compatible with each other (e.g., you don't have 2 instatoms in the same group that each make
            statements about a Dose, or each make a statement about Site. It is ok if they all use the same Directive and the same Anaphora). And each group
            has to provide a maximal cover of strucs in the struc_list, so that we can't add any more compatible atoms to the group. We do this via a recursive procedure
            find_maximal_instatom_groups_recursively.
        6.  Change these maximal groups to Versions and return the sorted list of Versions.

        """

        parse = instruction.parse
        strucs_labels_string = get_struc_labels(struc_list, delimiter = '|', omit_spaces_and_punctuation = True)
        if not 'DIRECTIVE' in strucs_labels_string and instruction.primary_directive:
            struc_list = [instruction.primary_directive] + struc_list
            strucs_labels_string += ' DIRECTIVE'
        struc_list = reorder_strucs_in_canonical_order_and_omit_punctuation(struc_list)
        compatible_atoms = get_compatible_atoms(struc_list)

        # Now find an "anaphorized" version of the struc_list where Dose (if it exists) is replaced with Anaphora "this medicine"
        if 'DIRECTIVE' in strucs_labels_string and 'ANAPHORA' not in strucs_labels_string:
            anaphorized_struc_list = struc_list[:]
            if 'DOSE' in strucs_labels_string:
                for struc in anaphorized_struc_list:
                    if struc.label == 'DOSE':
                        break
                anaphorized_struc_list.remove(struc)
            for i, struc in enumerate(anaphorized_struc_list):
                if struc.label == 'DIRECTIVE':
                    break
            if 'SPECIFIC_DAY' in strucs_labels_string:
                anaphora = Anaphora('this_dose', [])
            else:
                anaphora = Anaphora('this_medicine', [])
            anaphorized_struc_list.insert(i+1, anaphora)
        elif 'ANAPHORA' in strucs_labels_string:
            anaphorized_struc_list = struc_list
        else:
            anaphorized_struc_list = []

        #print('\n---------Event struc itself: \Keys %s\nAnaphorized Keys:%s\n-------\n' % (' '.join([str(get_key(struc)) for struc in struc_list]), ' '.join([str(get_key(struc)) for struc in anaphorized_struc_list])))
        matched_atoms = []
        # Match each compatible atom to the corresponding sub-sequence of the  sig strucs and see if the quality of the match is good enough.
        for atom in compatible_atoms:
            atom_strucs = reorder_strucs_in_canonical_order_and_omit_punctuation(atom.parses[0].strucs)
            atom_strucs_labels = [struc.label for struc in atom_strucs]
            #print('%d  Raw atom: %s  Labels: %s' % (compatible_atoms.index(atom), atom.raw_sig, ' '.join(atom_strucs_labels)))
            if 'ANAPHORA' not in atom_strucs_labels:
                # Try to match the regular strucs list
                struc_list_to_match = struc_list
            else:
                # Now try to match the anaphorized_struc_list to the atom
                struc_list_to_match = anaphorized_struc_list
            matching_strucs_lists = find_portion_of_struc_list_corresponding_to_atom(struc_list_to_match, atom_strucs, atom_strucs_labels)
            for matching_strucs in matching_strucs_lists:
                match_result = match_struc_list_to_atom(matching_strucs, atom)
                if match_result:
                    (quality_match_score, numeric_mappings) = match_result
                    instatom = InstAtom(atom, matching_strucs, quality_match_score, numeric_mappings)
                    matched_atoms.append(instatom)

        maximal_versions_list = group_matching_instatoms_into_versions(parse, matched_atoms)

        return maximal_versions_list

    def group_matching_instatoms_into_versions(parse, matched_instatoms_list):
        """ Takes a list of InstAtoms that all could be used for translation of portions of the sig, and returns the sorted list of Versions
            that comprize the maximal groups of compatible instatoms.

            Algorithm:
        1.  Prune matched_instatoms_list so that if two atoms from matched_instatoms_list match the same set of strucs, we keep only the one that matches best. We do this to keep the
            computational costs of the next step (the recursion to find maximal groups of these matched instatoms that are compatible with each other and
            provide maximal cover of strucs in the struc_list). The list of these instatoms is best_matched_atoms.
        2.  Group these matching Instatoms from best_matched_atoms into maximal groups maximal_instatom_groups. A maximal groups of these matched instatoms
            is a list of instatoms so that any two of them are compatible with each other (e.g., you don't have 2 instatoms in the same group that each make
            statements about a Dose, or each make a statement about Site. It is ok if they all use the same Directive and the same Anaphora). And each group
            has to provide a maximal cover of strucs in the struc_list, so that we can't add any more compatible atoms to the group. We do this via a recursive procedure
            find_maximal_instatom_groups_recursively.
        3.  Change these maximal groups to Versions and return the sorted list of Versions.
        """

        # For each matched_instatoms_list list (a Frame) from the sig, select the best quality atom that matches it.
        matched_instatoms_list.sort(key = lambda instatom: instatom.quality_match_score, reverse = True)
        matched_instatoms_list.sort(key = lambda instatom: instatom.atom.parses[0].flatten())
        best_matched_instatoms = []
        strucs_already_matched = []
        for instatom in matched_instatoms_list:
            if instatom.matching_strucs not in strucs_already_matched:
                strucs_already_matched.append(instatom.matching_strucs)
                best_matched_instatoms.append(instatom)

        # Now find maximal groups of InstAtoms to form versions.
        # Sorting instatoms from worst rank to best seems to lower the number of loops in find_maximal_instatom_groups_recursively()
        best_matched_instatoms.sort(key = lambda instatom: instatom.get_raw_rank())
        for rank, instatom in enumerate(best_matched_instatoms):
            instatom.rank = rank

        maximal_instatom_groups = find_maximal_instatom_groups_recursively(best_matched_instatoms)
        parse.debug_max_recursive_calls = max(parse.debug_max_recursive_calls, debug_recursive_calls[0])
        if debug:
            #print('\n\n-----Number of recursive calls to find_maximal_instatom_groups_recursively() is %d\n' % debug_recursive_calls[0])
            pass
        debug_recursive_calls[0] = 0

        empty_version = Version(parse)
        maximal_versions = []
        for instatom_group in maximal_instatom_groups:
            instatom_group.sort(key = lambda instatom: instatom.rank, reverse = True)
            version = empty_version.add_instatoms(instatom_group)
            maximal_versions.append(version)

        maximal_versions.sort(key = lambda version: version.quality_score(), reverse = True)

        if debug:
            for i, version in enumerate(maximal_versions):
                instatoms = version.instatoms
                txt = '|'.join([instatom.atom.parses[0].flatten() for instatom in instatoms])
                #print('--------- VERSION %d of %d: %s' % (i, len(maximal_versions), txt))

        return maximal_versions



    def find_maximal_instatom_groups_recursively(instatoms):
        """ Given a list of InstAtoms which match certain core strucs of the sig, find the maximal subsets of these that don't intersect on these core strucs.
            Return the list of these maximal subsets (actually, a list of maximal lists of instatoms).

            A maximal subset is one which can't be extended by any other instatom in the instatoms list.
            No two instatoms in one version (maximal subset) can have any core strucs from the sig in common because you can't state 2 different things about the same
            struc. The only exception (which are removed from the core strucs) are Directive and Anaphora, which can be repeated in all instatoms. So every 2 instatoms
            from the same max version have non-intersecting core_matched_strucs.

            Algorithm:
            -   Pick a Head Instatom (by default, the first one in whatever sort order provided in input).
            -   Find all instatoms in the universe incompatible with the Head (i.e. they have a core struc in common). Call this list incompatible_with_head. That list includes Head itself.
                -   Claim: Every maximal subset of instatoms has to have an element in incompatible_with_head (if a group doesn't you could add head to it, which means it was not maximal).
            -   For each incompatible_instatom in incompatible_with_head, find all the elements compatible with incompatible_instatom. Recurse on that list (minus preceding
                elements of incompatible_with_head) as the new universe, and append the element incompatible_instatom to each of the resultant maximal groups.
            -   This guarantees that we found all maximal groups and counted them exactly once. Proof: Consider a maximal group M. It has to have an element in incompatible_with_head,
                otherwise we could expand the group by adding head to it. Let X be the first element in incompatible_with_head that is in M.
                Then all elements of M - X are incompatible with X, and M contains none of the predecessors of X in incompatible_with_head, so M - X was included in the main step.
                And it was not included in any subsequent steps because they would have omitted X from the universe.

        """

        if not instatoms:
            return []
        elif len(instatoms) == 1:
            # If there is only 1 instatom, then there is only one maximal sublist: the list of length 1.
            return [instatoms]

        debug_recursive_calls[0] += 1

        maximal_instatom_groups = []

        head = instatoms[0]
        universe = instatoms[:]

        # compatibles are instatoms that don't haven any core_matched_strucs in common, i.e. the pairwise intersection of self.core_matched_strucs is empty.
        incompatible_with_head = [instatom for instatom in universe if head.is_incompatible_with(instatom)] # incompatible_with_head includes the head

        for incompatible_instatom in incompatible_with_head:
            universe.remove(incompatible_instatom)
            compatibles = incompatible_instatom.get_compatible_instatoms(universe)
            compatible_instatom_groups = find_maximal_instatom_groups_recursively(compatibles)
            if compatible_instatom_groups:
                new_groups = [[incompatible_instatom] + instatom_group for instatom_group in compatible_instatom_groups]
            else:
                new_groups = [[incompatible_instatom]]

            maximal_instatom_groups += new_groups

        return maximal_instatom_groups



    def process_event(event, is_initial_event, omit_these_strucs, include_then_chrono_in_event, omit_and_conj_in_event, process_schedule_separately):

        schedule = event.schedule
        instruction = schedule.instruction
        parse = instruction.parse
        empty_version = Version(parse)

        event_components = event.get_recursive_componenets()
        schedule_components = schedule.get_prop_strucs()

        and_conj_struc = None
        then_chrono_struc = ThenChrono([]) if include_then_chrono_in_event else None

        if is_initial_event:
            if process_schedule_separately:
                all_components = event_components
            else:
                all_components = event_components + schedule_components
            if include_then_chrono_in_event:
                all_components.insert(0, then_chrono_struc)
        else:
            if omit_and_conj_in_event:
                all_components = event_components
            else:
                and_conj_struc = AndConj([])
                all_components = [and_conj_struc] + event_components

        if omit_these_strucs:
            all_components = [struc for struc in all_components if struc not in omit_these_strucs]

        if not all_components:
            return []

        new_versions = []

        # Main Process.
        versions = match_sem_chunk_to_dictionary(all_components, instruction)
        new_versions += versions

        # If there are uncovered strucs in the best scoring versions, try to change some things, e.g. adjust DIRECTIVE if it seems to conflict with FORM or Vehicle
        uncovered_label2struc_list = find_strucs_not_covered_by_top_versions(all_components, versions, max_number_of_versions = 2)
        versions = alter_directive_if_needed(all_components, event, versions, uncovered_label2struc_list)
        new_versions += versions
        if 'TIMING' in uncovered_label2struc_list:
            versions = case_timing_not_covered(all_components, event, new_versions, uncovered_label2struc_list)
            new_versions += versions
        if 'ROUTE' in uncovered_label2struc_list:
            versions = case_route_not_covered(all_components, event, new_versions, uncovered_label2struc_list)
            new_versions += versions


        # penalize versions for absence of Then_Chrono and/or And_Conj
        # if is a non-initial schedule, then the initial event should have 'THEN_CHRONO'.
        # if is a non-initial event, then should have 'AND_CONJ'
        for version in new_versions:
            version.assess_quality_penalty_for_event(event = event, is_initial_event = is_initial_event)
            if and_conj_struc and and_conj_struc not in version.matched_strucs:
                version.quality_penalty += 3
            if then_chrono_struc and then_chrono_struc not in version.matched_strucs:
                version.quality_penalty += 3

        new_unique_versions = Version.prune_and_sort_versions(new_versions)
        return new_unique_versions

    def find_strucs_not_covered_by_top_versions(struc_list, versions, max_number_of_versions = None):
        """ Finds all strucs in struc_list that are not covered by any of the top versions in the versions list. Returns a dictionary
            label2struc_list which maps labels of such uncovered versions to the list of version strucs.

            Because struc_list sometimes has different (alternative) strucs to those in the Sem, we can't use version.strucs_not_covered() function
            but have to see if each struc in struc_list is in any version.matched_strucs
        """


        if max_number_of_versions is None:
            max_number_of_versions = 2

        label2struc_list = defaultdict(list)
        if not versions:
            uncovered_strucs = set(struc_list)
        else:
            versions.sort(key = lambda version: version.quality_score(), reverse = True)
            covered_strucs = set()
            for version in versions[:max_number_of_versions]:
                covered_strucs |= version.matched_strucs
            uncovered_strucs = set(struc_list) - covered_strucs

        for struc in uncovered_strucs:
                label2struc_list[struc.label].append(struc)

        return label2struc_list

    def case_timing_not_covered(all_components, event, versions, uncovered_label2struc_list):
        """ If some Timing is not covered by any top version, tries to alter Timing.

        Deals with the cases:
        -   "nightly" (e.g. "take 1 tablet by mouth nightly" and tries to change "night" to "bedtime" if "night" is not covered by any version.
        -   offset in timing, e.g. "take 1/2 hour before breakfast", but the dicitonary may have "30 minutes before breakfast". So if timing is not matched, we try to change
            hour to minutes and vice versa
        -   offset in timing is not vital and not in dictionary, e.g. "take within 30 mins of onset of headache" is not in the dictionary.

        """

        if 'TIMING' not in uncovered_label2struc_list:
            return []

        uncovered_timings = uncovered_label2struc_list['TIMING']
        new_versions = []
        night_timing = [timing for timing in uncovered_timings if timing.landmark == 'night']
        offset_timing = [timing for timing in uncovered_timings if timing.offset]
        if night_timing:
            night_timing = night_timing[0]
            night_timing.landmark = 'bedtime'
            new_versions = match_sem_chunk_to_dictionary(all_components, event.schedule.instruction)
            if new_versions:
                # Check if the best new version now contains TIMING.
                new_versions.sort(key = lambda version: version.quality_score(), reverse = True)
                new_uncovered_label2struc_list = find_strucs_not_covered_by_top_versions(all_components, new_versions, max_number_of_versions = 1)
            if not new_versions or 'TIMING' in new_uncovered_label2struc_list:
                night_timing.landmark = 'night' # revert to old Timing.
                new_versions = []
        elif offset_timing:
            offset_timing = offset_timing[0]
            old_quant = offset_timing.offset.quant
            old_time_unit_value = offset_timing.offset.time_unit.value
            if old_time_unit_value == 'minute':
                new_time_unit_value = 'hour'
                if old_quant.num_type == 'range':
                    old_quant_value = old_quant.low.value
                elif old_quant.num_type == 'var':
                    return []
                else:
                    old_quant_value = old_quant.value
                new_quant_value = old_quant_value / 60.0
                if new_quant_value == int(new_quant_value):    # integral number of hours
                    new_quant = Quant('int', int(new_quant_value))
                elif old_quant_value == 5:
                    new_quant = FracQuant(0, 1, 12)
                elif old_quant_value == 10:
                    new_quant = FracQuant(0, 1, 6)
                elif old_quant_value == 15:
                    new_quant = FracQuant(0, 1, 4)
                elif old_quant_value == 20:
                    new_quant = FracQuant(0, 1, 3)
                elif old_quant_value == 30:
                    new_quant = FracQuant(0, 1, 2)
                elif old_quant_value == 40:
                    new_quant = FracQuant(0, 2, 3)
                elif old_quant_value == 45:
                    new_quant = FracQuant(0, 3, 4)
                elif old_quant_value == 90:
                    new_quant = Quant('decimal', 1.5)
                else:
                    new_quant = Quant('decimal', round(new_quant_value, 2))
            elif old_time_unit_value == 'hour':
                new_time_unit_value = 'minute'
                if old_quant.num_type == 'range':
                    old_quant_value = old_quant.low.value
                elif old_quant.num_type == 'var':
                    return []
                else:
                    old_quant_value = old_quant.value
                new_quant = Quant('int', int(old_quant_value * 60))
            else:
                return []
            offset_timing.offset.quant = new_quant
            offset_timing.offset.time_unit.value = new_time_unit_value
            new_versions = match_sem_chunk_to_dictionary(all_components, event.schedule.instruction)
            if new_versions:
                # Check if the best new version now contains TIMING.
                new_versions.sort(key = lambda version: version.quality_score(), reverse = True)
                new_uncovered_label2struc_list = find_strucs_not_covered_by_top_versions(all_components, new_versions, max_number_of_versions = 1)
            if not new_versions or 'TIMING' in new_uncovered_label2struc_list:
                offset_timing.offset.quant =  old_quant     # revert to old Timing.
                offset_timing.offset.time_unit.value = old_time_unit_value
                new_versions = []
                if offset_timing.landmark == 'onset_of_headache':
                    old_time_unit = offset_timing.offset.time_unit
                    offset_timing.offset = None
                    new_versions = match_sem_chunk_to_dictionary(all_components, event.schedule.instruction)
                    if new_versions:
                        new_versions.sort(key = lambda version: version.quality_score(), reverse = True)
                        new_uncovered_label2struc_list = find_strucs_not_covered_by_top_versions(all_components, new_versions, max_number_of_versions = 1)
                    if not new_versions or 'TIMING' in new_uncovered_label2struc_list:
                        offset_timing.offset = TimeInterval(quant = old_quant, time_unit = old_time_unit) # revert to old Timing.
                        new_versions = []

        return new_versions

    def case_route_not_covered(all_components, event, old_versions, uncovered_label2struc_list):
        """ If some Route is not covered by any top version, tries to force-add a simple atom to match it.

        Example:
        "Topically" is usually used with "apply" directive, but sometimes there is "use topically", which we don't have in the dictionary, so we force-match it to "apply this medicine topically"
        However, we don't want to change the directive in the event or in the instruction.primary_directive because "apply" matches the other things (such as DOSE) just fine.
        """

        if 'ROUTE' not in uncovered_label2struc_list:
            return []

        uncovered_route = uncovered_label2struc_list['ROUTE'][0]
        current_directive = event.schedule.instruction.primary_directive
        expected_directive_value = Route.value_to_likely_directive.get(uncovered_route.value, '')
        new_versions = []
        if current_directive and current_directive.value != expected_directive_value:
            expected_directive = Directive(expected_directive_value)
            anaphora = Anaphora(value = 'this_medicine', constituents = [])
            struc_list = [expected_directive, uncovered_route]
            candidate_atoms = get_structurally_similar_list_of_atoms_in_dict(struc_list)    # Try "Apply topically")
            new_versions = try_adding_struc_list_to_versions_list(struc_list, old_versions, restrict_to_these_atoms_list = candidate_atoms)

            if not new_versions:    # If without Anaphora didn't work, try with Anaphora ("Apply this medicine topically")
                struc_list = [expected_directive, anaphora, uncovered_route]
                candidate_atoms = get_structurally_similar_list_of_atoms_in_dict(struc_list)
                new_versions = try_adding_struc_list_to_versions_list(struc_list, old_versions, restrict_to_these_atoms_list = candidate_atoms)
        return new_versions



    def alter_directive_if_needed(struc_list, event, old_versions, uncovered_label2struc_list):
        """ If the top 2 versions didn't cover some top Sems (DOSE, Site, Vehicle), returns additional versions by trying special cases. Returns list of new_versions, if any.

        If we didn't cover dose, it could be that the directive is unexpected.
        Example: "USE 1 SUPPOSOTORY RECTALLY AS DIRECTED" -- is "USE" but all atoms for form = suppository have "Insert suppository". So try the alternative directive.

        Similarly with VEHICLE (nebulizer goes with inhale in dictionary but is used often with "use").
        Similarly for SITE ("eye" and "ear" go with instill, not use or place).
        """

        instruction = event.schedule.instruction
        directive_in_event = [struc for struc in struc_list if struc.label == 'DIRECTIVE']
        if directive_in_event:
            directive_in_event = directive_in_event[0]
            directive = directive_in_event
        else:
            directive = instruction.primary_directive
        old_directive_value = directive.value if directive else ''

        all_new_versions = []

        for label in ('DOSE', 'VEHICLE', 'SITE'):
            if label in uncovered_label2struc_list:
                labeled_struc = uncovered_label2struc_list[label][0]
            else:
                continue

            if label == 'DOSE':
                expected_directive_value = Form.value_to_likely_directive.get(labeled_struc.form.value, '')
            elif label == 'VEHICLE':
                expected_directive_value = Vehicle.value_to_likely_directive.get(labeled_struc.value, '')
            elif label == 'SITE':
                expected_directive_value = Site.value_to_likely_directive.get(labeled_struc.value, '')

            if directive.value == 'swish and spit':
                expected_directive_value = 'rinse'

            if expected_directive_value and old_directive_value != expected_directive_value:
                directive.value = expected_directive_value
                old_primary_directive_value = instruction.primary_directive.value if instruction.primary_directive else ''
                new_versions = match_sem_chunk_to_dictionary(struc_list, instruction)
                if new_versions:
                    # Check if the best new version now contains the label.
                    new_uncovered_label2struc_list = find_strucs_not_covered_by_top_versions(struc_list, new_versions, max_number_of_versions = 1)
                    if label not in new_uncovered_label2struc_list: # The top new version covers this struc
                        instruction.primary_directive = directive
                        directive.rules_used.append('modified_in_alter_directive_if_needed_for_dictionary_matching_0')
                        all_new_versions += new_versions
                if not new_versions or label in new_uncovered_label2struc_list:    # Nothing found with this attempted change. Undo change of directive in the Sem
                    directive.value = old_directive_value

            if label == 'DOSE' and labeled_struc.form.value == 'inhalation':
                # Sigs such as "use 2 inhalations" need to have both Form changed (to puffs) and directive to "inhale"
                directive.value = 'inhale'
                labeled_struc.form.value = 'puff'
                new_versions = match_sem_chunk_to_dictionary(struc_list, instruction)
                if new_versions:
                    # Check if the best new version now contains DOSE.
                    new_uncovered_label2struc_list = find_strucs_not_covered_by_top_versions(struc_list, new_versions, max_number_of_versions = 1)
                    if label not in new_uncovered_label2struc_list: # The top new version covers this struc
                        instruction.primary_directive = directive
                        instruction.primary_form = labeled_struc.form
                        directive.rules_used.append('modified_in_alter_directive_if_needed_for_dictionary_matching_1')
                        labeled_struc.form.rules_used.append('modified_in_alter_directive_if_needed_for_dictionary_matching_1')
                        all_new_versions += new_versions
                if not new_versions or label in new_uncovered_label2struc_list:         # Nothing found with this attempted change. Undo change of Directive and DOSE in the Sem
                    directive.value = old_directive_value
                    labeled_struc.form.value = 'puff'

            elif label == 'DOSE' and labeled_struc.form.value == 'application' and event.site and event.site.value in ('eye', 'ear'):
                # Sigs such as "1 application in both eyes" need to have both Form changed (to drops) and directive to "instill"
                directive.value = 'instill'
                labeled_struc.form.value = 'drop'
                new_versions = match_sem_chunk_to_dictionary(struc_list, instruction)
                if new_versions:
                    # Check if the best new version now contains DOSE.
                    new_uncovered_label2struc_list = find_strucs_not_covered_by_top_versions(struc_list, new_versions, max_number_of_versions = 1)
                    if label not in new_uncovered_label2struc_list: # The top new version covers this struc
                        instruction.primary_directive = directive
                        instruction.primary_form = labeled_struc.form
                        directive.rules_used.append('modified_in_alter_directive_if_needed_for_dictionary_matching_2')
                        labeled_struc.form.rules_used.append('modified_in_alter_directive_if_needed_for_dictionary_matching_2')
                        all_new_versions += new_versions
                if not new_versions or label in new_uncovered_label2struc_list:       # Nothing found with this attempted change. Undo change of Directive and DOSE in the Sem
                    directive.value = old_directive_value
                    labeled_struc.form.value = 'application'

            elif label == 'DOSE' and labeled_struc.form.value == 'gram' and event.substrate:
                # Sigs such as "dissolve 17 grams of powder in 8 ounces of water" need to have directive changed to "mix"
                directive.value = 'mix'
                new_versions = match_sem_chunk_to_dictionary(struc_list, instruction)
                if new_versions:
                    # Check if the best new version now contains DOSE.
                    new_uncovered_label2struc_list = find_strucs_not_covered_by_top_versions(struc_list, new_versions, max_number_of_versions = 1)
                    if label not in new_uncovered_label2struc_list: # The top new version covers this struc
                        instruction.primary_directive = directive
                        directive.rules_used.append('modified_in_alter_directive_if_needed_for_dictionary_matching_3')
                        all_new_versions += new_versions
                if not new_versions or label in new_uncovered_label2struc_list:       # Nothing found with this attempted change. Undo change of Directive and DOSE in the Sem
                    directive.value = old_directive_value

            #elif label == 'SITE' and labeled_struc.value == 'eyelid':
            #    # Sigs such as "apply in both eyelids", if "eyelid" is not in dictionary, "eye" may be. So need to have both Site changed (to "eye") and directive to "instill"
            #    directive.value = 'instill'
            #    site = labeled_struc
            #    site.value = 'eye'
            #    if event.dose:
            #        old_form_value = event.dose.form.value
            #        event.dose.form.value = 'drop'
            #    new_versions = match_sem_chunk_to_dictionary(struc_list, instruction)
            #    if new_versions:
            #        # Check if the best new version now contains SITE.
            #        new_uncovered_label2struc_list = find_strucs_not_covered_by_top_versions(struc_list, new_versions, max_number_of_versions = 1)
            #        if label not in new_uncovered_label2struc_list: # The top new version covers this struc
            #            instruction.primary_directive = directive
            #            instruction.primary_form = event.dose.form if event.dose else Form('drop', plurality  = 'plurality_either', constituents = [])
            #            directive.rules_used.append('modified_in_alter_directive_if_needed_for_dictionary_matching_4')
            #            instruction.primary_form.rules_used.append('modified_in_alter_directive_if_needed_for_dictionary_matching_4')
            #            all_new_versions += new_versions
            #    if not new_versions or label in new_uncovered_label2struc_list:       # Nothing found with this attempted change. Undo change of Directive, SITE, and DOSE in the Sem
            #        directive.value = old_directive_value
            #        site.value = 'eyelid'
            #        if event.dose:
            #            event.dose.form.value = old_form_value

            elif label == 'SITE' and labeled_struc.value == 'nostril' and event.dose and event.dose.form.value == 'puff':
                # Sigs such as "spray two puffs to each nostril", where the dictionary uses form = "spray" (e.g. "use 2 sprays to each nostril" or "spray 2 times to each nostril").
                # So need to have both Dose.form changed (to "spray") and directive to "spray" or "use"
                event.dose.form.value = 'spray'
                directive.value = 'spray'
                new_versions = match_sem_chunk_to_dictionary(struc_list, instruction)
                if new_versions:
                    # Check if the best new version now contains SITE.
                    new_uncovered_label2struc_list = find_strucs_not_covered_by_top_versions(struc_list, new_versions, max_number_of_versions = 1)
                    if label not in new_uncovered_label2struc_list: # The top new version covers this struc
                        instruction.primary_directive = directive
                        instruction.primary_form = event.dose.form
                        directive.rules_used.append('modified_in_alter_directive_if_needed_for_dictionary_matching_5')
                        event.dose.form.rules_used.append('modified_in_alter_directive_if_needed_for_dictionary_matching_5')
                        all_new_versions += new_versions
                    else:   # Try "use" instead of "spray" for directive
                        directive.value = 'use'
                        new_versions = match_sem_chunk_to_dictionary(struc_list, instruction)
                    if new_versions:
                        # Check if the best new version now contains SITE.
                        new_uncovered_label2struc_list = find_strucs_not_covered_by_top_versions(struc_list, new_versions, max_number_of_versions = 1)
                        if label not in new_uncovered_label2struc_list: # The top new version covers this struc
                            instruction.primary_directive = directive
                            instruction.primary_form = event.dose.form
                            directive.rules_used.append('modified_in_alter_directive_if_needed_for_dictionary_matching_6')
                            event.dose.form.rules_used.append('modified_in_alter_directive_if_needed_for_dictionary_matching_6')
                            all_new_versions += new_versions
                    if not new_versions or label in new_uncovered_label2struc_list:       # Nothing found with this attempted change. Undo change of Directive, DOSE in the Sem
                        directive.value = old_directive_value
                        event.dose.form.value = 'puff'


        return all_new_versions



    def process_schedule_proper(schedule, insert_then_chrono, exclude_duration_from_schedule):
        """ Creates a list of versions that transduce the strucs that are semantic components of the Schedule proper (i.e. not of any of the Events that
            are part of the schedule. This is used only if can't process the schedule proper components with the first event of the schedule (e.g.
            because we have multiple events in the schedule, and we can't talk about Duration before talking about the second event (e.g.
            instead of "Take 2 tabs in the AM for 3 days. Also: take 1 tab at bedtime" you want to say "Take 2 tabs in the AM. Also: take 1 tab at bedtime. Take this medicine for 3 days."

        """

        instruction = schedule.instruction
        schedule_components = schedule.get_prop_strucs()
        if insert_then_chrono:
            then_chrono_struc =  ThenChrono([])
            schedule_components.insert(0, then_chrono_struc)
        if exclude_duration_from_schedule and schedule.duration:
            schedule_components.remove(schedule.duration)
        versions = match_sem_chunk_to_dictionary(schedule_components, instruction)
        # If the only thing the versions matched was Then_Chrono, they should be deleted
        pruned_versions = []
        for version in versions:
            if insert_then_chrono and len(version.matched_strucs) == 1 and then_chrono_struc in version.matched_strucs:
                # The version transduces the whole schedule into just "After that:". Remove this version.
                continue
            else:
                pruned_versions.append(version)
        unique_versions = Version.prune_and_sort_versions(pruned_versions)
        return unique_versions

    def process_schedule(schedule, is_initial_schedule):
        """ Creates a list of brand new versions that transduce this schedule. Does not try to process strucs that belong to the Instruction level --
            they are processed separately.

            If not initial schedule, don't process duration or calendar_event, because they were already appended before this.


        Contraint:
        A case such as this:        TAKE ONE TABLET BY MOUTH TWICE DAILY BEFORE BREAKFAST AND DINNER
        should have at least 2 versions:
            Ver 1: Take 1 tablet by mouth before breakfast. Also: take this medicine before dinner. Take this medicine 2 times a day. Take this medicine every day.
            Ver 2: Take 1 tablet by mouth 2 times a day. Take this medicine every day. Take this medicine before breakfast. Also: take this medicien before dinner.
            Ver 3: TAKE 1 TABLET BY MOUTH 2 TIMES A DAY. TAKE THIS MEDICINE BEFORE BREAKFAST AND BEFORE DINNER.
            We should NOT have "Take 1 tablet by mouth 2 times a day before breakfast. Also take this medicine before dinner" because this falsely implies that Freq pertains to "before breakfast."
        """

        empty_version = Version(parse)

        # If there are no events in the schedule (e.g. this is a second Schedule and it only has Freq ("Take twice daily for the first week then once a day")
        # Then we need to signal to process schedule components later.
        process_schedule_separately = False if schedule.events else True
        include_then_chrono_in_event = False
        exclude_duration_from_schedule = False

        new_versions = []
        for event_num, event in enumerate(schedule.events):
            if event_num == 0:
                # If we have multiple events, then schedule components such as DURATION, FREQ should go after all events are processed.
                # E.g. you don't want to have "Take 2 tabs in the AM for 3 days. Also: take 1 tab at bedtime".

                if len(schedule.events) <= 1:
                    process_schedule_separately = False
                elif len(schedule.events) == 2 and not schedule.events[0].timing and schedule.events[1].timing and (schedule.freq or schedule.periodicity):
                    # An exception is a schedule that has Freq or Periodicity but no Timing in the first event and only Timing in the second event.
                    # Then the Freq of the Schedule really pertains to the First Event and should be processed with it.
                    # E.g. "apply rectally 2-3 times a day and after each bowel movement". Here "2-3 times a day" is part of the first Event, and "After bowel movement" is a second Event.
                    # Another example: "apply one patch one time daily (remove at bedtime)"
                    process_schedule_separately = False
                else:
                    process_schedule_separately = True

                # For initial_event in a non-initial schedule, we need to preface it with "Then". Prima facie, the easiest thing to do
                # is put "Then:" before processing the Event. But that's not always optimal. Sometimes
                # this "Then_chrono" is accomodated within a full-bodied atom, e.g. THEN TAKE 1 TEASPOON ON DAYS <<NUM_0>> - <<NUM_1>>.
                # so we want to simply signal to () to include Then in the list of strucs that have to be matched by InstAtoms.
                # But if there is a Duration struc in this non-initial Schedule, we may want to preface the event with "Then Duration:"
                # (e.g., "THEN FOR THE NEXT <<NUM_0>> DAYS:"). So in these cases we should try both. In the latter case (when we pull out the duration upfront)
                # we need to make sure that we yank it out when we process the schedule strucs proper in process_schedule_proper()

                if is_initial_schedule:
                    include_then_chrono_in_event = False
                else:
                    include_then_chrono_in_event = True

                if not is_initial_schedule and process_schedule_separately and schedule.duration:
                    original_offset = schedule.duration.offset
                    if not original_offset:
                        schedule.duration.offset = 'next'
                    then_duration = [ThenChrono([]), schedule.duration, Struc(label = ':')]
                    candidate_atoms = get_structurally_similar_list_of_atoms_in_dict(then_duration)
                    then_duration_versions = try_adding_struc_list_to_versions_list(then_duration, [empty_version], restrict_to_these_atoms_list = candidate_atoms)
                    if then_duration_versions:
                        include_then_chrono_in_event = False
                        exclude_duration_from_schedule = True
                    else:
                        schedule.duration.offset = original_offset

                event_0_versions = process_event(event, is_initial_event = True,
                                                 omit_these_strucs = set(),
                                                 include_then_chrono_in_event = include_then_chrono_in_event,
                                                 omit_and_conj_in_event = False,
                                                 process_schedule_separately = process_schedule_separately)
                if not event_0_versions:
                    # If we can't process event_0, don't process the other events by sticking "Also:" in front of them.
                    break
                new_versions = event_0_versions
            else:   # post-initial events
                new_event_versions = []

                if event.directive and event.directive.value == 'remove' and event.timing:
                    # Special case for Remove directive, because it makes no sense to say "Also: remove in the evening".
                    # "Remove" is not an additional event but a reversal of the first Event.
                    omit_and_conj_in_event = True
                elif event.directive and event.directive.value in ('drink', 'inject') and schedule.events[0].directive and schedule.events[0].directive.value in ('mix', 'dissolve'):
                    # This second event is not a second thing to do during the day but is part of the fulfillment of the first event, so don't add "Also:"
                    omit_and_conj_in_event = True
                else:
                    omit_and_conj_in_event = False

                versions = process_event(event, is_initial_event = False,
                                         omit_these_strucs = set(),
                                         include_then_chrono_in_event = False,
                                         omit_and_conj_in_event = omit_and_conj_in_event,
                                         process_schedule_separately = True)
                versions = Version.concatenate_versions(new_versions, versions)
                new_versions += versions

        if process_schedule_separately:
            if (is_initial_schedule or schedule.events) and not include_then_chrono_in_event:
                # We only need to insert "Then" into the Schedule verion if this a non-initial schedule which has no events (e.g. just Freq or Duration)
                insert_then_chrono = False
            else:
                insert_then_chrono = True
            schedule_versions = process_schedule_proper(schedule, insert_then_chrono = insert_then_chrono, exclude_duration_from_schedule = exclude_duration_from_schedule)
            new_versions = Version.concatenate_versions(new_versions, schedule_versions)
            if exclude_duration_from_schedule:
                # If we pulled out "Then Duration" (e.g. "then for the next 4 days:") from Schedule, include these upfront.
                new_versions = Version.concatenate_versions(then_duration_versions, new_versions)

        # Special cases
        if len(schedule.events) == 2:
            event_0 = schedule.events[0]
            event_1 = schedule.events[1]
            if event_0.timing and event_1.timing and len(event_1.get_recursive_componenets()) == 1:
                # Special case to deal with sigs that have 2 timing components with the same Dose: because in the dictionary we have serveral instances
                # of "DIRECTIVE ANAPHORA TIMING AND_CONJ TIMING" and "DIRECTIVE DOSE TIMING AND_CONJ TIMING" atoms, we try to see if they can be matched.
                directive = instruction.primary_directive
                anaphora = Anaphora(value = 'this_medicine', constituents = [])
                and_conj = AndConj([])
                if event_0.dose:
                    # Try "DIRECTIVE DOSE TIMING AND_CONJ TIMING"
                    dose = event_0.dose
                    struc_list = [directive, dose, event_0.timing[0], and_conj, event_1.timing[0]]
                    versions = process_special_2event_timing_case(event_0, struc_list, special_event_strucs_set = set(event_0.timing + [dose]),
                                                                  is_initial_schedule = is_initial_schedule,
                                                                  process_schedule_separately = True)   # We process_schedule_separately (meaning after both events)
                                                                                                        # because we are covering the dose here, so struc_list has to go upfront in the final version.
                    new_versions += versions
                # Try "DIRECTIVE ANAPHORA TIMING AND_CONJ TIMING"
                struc_list = [directive, anaphora, event_0.timing[0], and_conj, event_1.timing[0]]
                versions = process_special_2event_timing_case(event_0, struc_list, special_event_strucs_set = set(event_0.timing),
                                                              is_initial_schedule = is_initial_schedule,
                                                              process_schedule_separately = False)  # No need to process_schedule_separately because process_event() for event_0 will take care of all
                                                                                                    # the schedule components, and they can well go before this struc_list.
                new_versions += versions

        new_unique_versions = Version.prune_and_sort_versions(new_versions)
        return new_unique_versions

    def process_special_2event_timing_case(event_0, struc_list, special_event_strucs_set, is_initial_schedule, process_schedule_separately):

        empty_version = Version(parse)
        timing_versions = try_adding_struc_list_to_versions_list(struc_list, [empty_version])
        versions = []
        if timing_versions:
            event_0_versions = process_event(event_0,
                                             is_initial_event = True,
                                             include_then_chrono_in_event = not is_initial_schedule,
                                             omit_and_conj_in_event = False,
                                             omit_these_strucs = special_event_strucs_set,
                                             process_schedule_separately = process_schedule_separately)
            if process_schedule_separately:
                schedule_versions = process_schedule_proper(event_0.schedule, insert_then_chrono = False, exclude_duration_from_schedule = False)
            else:
                schedule_versions = []
            versions = Version.concatenate_versions(event_0_versions, timing_versions, schedule_versions)
        return versions

    def process_instruction_stucs_proper(instruction):
        """ Creates a list of brand new versions that transduce just the proper properties of the Instruction, i.e. Indication, As_Needed, See_Instructions, etc.
        """

        all_struc_list = instruction.get_prop_strucs()

        indications = instruction.indication
        as_needed = instruction.as_needed

        subsequent_indications_versions = []
        directive = instruction.primary_directive

        if indications and len(indications) > 1 and directive:
            # The only strucs that can repeat at the Instruction level are indications. The only possible matching atoms for indications
            # are DIRECTIVE ANAPHORA AS_NEEDED INDICATION, DIRECTIVE ANAPHORA INDICATION, DIRECTIVE AS_NEEDED INDICATION, DIRECTIVE INDICATION,
            # It is best stylistically to add "as_needed" to the first indication when there are several indications. So we deal with the first indication last
            empty_version = Version(parse)
            anaphora = Anaphora(value = 'this_medicine', constituents = [])
            for indication in indications[1:]:
                struc_list = [directive, indication]
                new_versions = try_adding_struc_list_to_versions_list(struc_list, [empty_version])
                struc_list = [directive, anaphora, indication]
                new_versions += try_adding_struc_list_to_versions_list(struc_list, [empty_version])
                if not new_versions and as_needed:
                    struc_list = [directive, as_needed, indication]
                    new_versions = try_adding_struc_list_to_versions_list(struc_list, [empty_version])
                    struc_list = [directive, anaphora, as_needed, indication]
                    new_versions += try_adding_struc_list_to_versions_list(struc_list, [empty_version])
                    if new_versions and as_needed in all_struc_list:
                        all_struc_list.remove(as_needed)
                if new_versions:
                    and_conj = [AndConj([])]
                    candidate_and_conj_atoms = get_structurally_similar_list_of_atoms_in_dict(and_conj)
                    and_conj_versions = try_adding_struc_list_to_versions_list(and_conj, [empty_version], restrict_to_these_atoms_list = candidate_and_conj_atoms)
                    subsequent_indications_versions = Version.concatenate_versions(subsequent_indications_versions, and_conj_versions, new_versions)
                all_struc_list.remove(indication)

        new_versions = match_sem_chunk_to_dictionary(all_struc_list, instruction)
        new_versions.sort(key = lambda version: version.quality_score(), reverse = True)
        # Check if we covered the indication in the best new_version (or if there are no versions found). If not, and as_needed is not in the sig, try adding it.
        if indications and (not new_versions or indications[0] in new_versions[0].strucs_not_covered()):
            if not as_needed:
                as_needed = AsNeeded([])
            additional_versions = match_sem_chunk_to_dictionary(all_struc_list + [as_needed], instruction)
            # But if the additional_versions still don't cover the indication, don't use them (they are likely just to contain "use as needed")
            additional_versions.sort(key = lambda version: version.quality_score(), reverse = True)
            if additional_versions and indications[0] in additional_versions[0].strucs_not_covered():
                pass
            else:
                new_versions += additional_versions
        all_new_versions = Version.concatenate_versions(new_versions, subsequent_indications_versions)
        return all_new_versions

    def process_special_duration_cases(instruction):
        """ Deals with special cases, where it extracts out of the Sem structure and processes one or multiple durations in custom ways.

            The expressions related to patches "apply 12 hours on then 12 hours off" are currently Sem'ed as two Schedules, one with Duration/on-12-hours,
            the other Duration/off-12-hours. But this is a non-intended use of "Duration" because the 12 hrs/on 12hrs/off cycle is happening within 1 day, and
            is repeated every day, so this is really one schedule.
        """

        def find_pair_on_off_duration_schedules(instruction):
            for sched_0_index in range(len(instruction.schedules) - 1):
                (schedule_0, schedule_1) = instruction.schedules[sched_0_index : sched_0_index + 2]
                duration_0 = schedule_0.duration
                duration_1 = schedule_1.duration
                if (duration_0 and duration_1 and duration_0.on_off == 'on' and duration_1.on_off == 'off' and
                    duration_0.time_interval.time_unit.value == duration_1.time_interval.time_unit.value and duration_0.time_interval.time_unit.value == 'hour'):
                    return (schedule_0, schedule_1)
            return None

        if instruction.schedules and len(instruction.schedules) > 1 and instruction.primary_directive.value == 'apply' and instruction.primary_form.value == 'patch':
            on_off_duration_schedules = find_pair_on_off_duration_schedules(instruction)
            if not on_off_duration_schedules:
                return []
            (schedule_0, schedule_1) = on_off_duration_schedules
            duration_0 = schedule_0.duration
            duration_1 = schedule_1.duration
            # we are dealing with "apply (patch) 12 hours on then 12 hours off"
            directive = instruction.primary_directive
            if schedule_1.events:
                directive_remove = schedule_1.events[0].directive
            else:   # should never happen, since schedule_1 should have 1 event with directive = remove
                return []
            then_chrono = ThenChrono([])
            struc_list = [directive, duration_0, then_chrono, directive_remove, duration_1]
            dose_included_event = None
            if schedule_0.events and schedule_0.events[0].dose:
                dose = schedule_0.events[0].dose
                struc_list.insert(1, dose)
                dose_included_event = schedule_0.events[0]
            elif schedule_1.events and schedule_1.events[0].dose:
                # case such as "apply 12 hours on and 12 hours off daily 1 patch"
                dose = schedule_1.events[0].dose
                struc_list.insert(1, dose)
                dose_included_event = schedule_1.events[0]
            elif instruction.schedules[0].events and  instruction.schedules[0].events[0].dose:
                # case where the on/off duration are not in the first schedule, e.g. "apply 1-3 patches topically once daily for 30 days 12 hours on and 12 hours off daily"
                dose = instruction.schedules[0].events[0].dose
                struc_list.insert(1, dose)
                dose_included_event = instruction.schedules[0].events[0]
            empty_version = Version(parse)

            candidate_atoms = get_structurally_similar_list_of_atoms_in_dict(struc_list)
            versions = try_adding_struc_list_to_versions_list(struc_list, [empty_version], restrict_to_these_atoms_list = candidate_atoms)
            if versions:
                # Check if the best new version now contains Duration. If it does, we then change the Sem by consolidating the 2 schedules into 1 and removing
                # both Duration
                uncovered_label2struc_list = find_strucs_not_covered_by_top_versions(struc_list, versions, max_number_of_versions = 1)
                if 'DURATION' not in uncovered_label2struc_list:
                    # The top new version covers Duration
                    schedule_1_properties = [prop for prop in schedule_1.substantive_properties if schedule_1.get_property_values(prop)]
                    for prop in schedule_1_properties:
                        # for every substantive component of schedule_1 proper, reassign it to schedule_0.
                        value = schedule_1.get_property_values(prop)
                        if type(value) == list:
                            for struc in value:
                                schedule_1.remove_property(prop, struc)
                                if not schedule_0.get_property_values(prop):
                                    schedule_0.add_property_value(prop, struc)
                                    struc.rules_used.append('*reassigned_sem_by_process_special_duration_cases()*')
                                else:
                                    struc.rules_used.append('*removed_from_sem_by_process_special_duration_cases()*')
                        else:
                            schedule_1.remove_property(prop, value)
                            if not schedule_0.get_property_values(prop):
                                schedule_0.add_property_value(prop, value)
                                value.rules_used.append('*reassigned_sem_by_process_special_duration_cases()*')
                            else:
                                value.rules_used.append('*removed_from_sem_by_process_special_duration_cases()*')
                    schedule_0.remove_property('duration', duration_0)
                    duration_0.rules_used.append('*removed_from_sem_by_process_special_duration_cases()*')
                    index_schedule_1 = instruction.schedules.index(schedule_1)
                    instruction.schedules[index_schedule_1:index_schedule_1 + 1] = []
                    if dose_included_event:
                        dose_included_event.dose.rules_used.append('*removed_from_sem_by_process_special_duration_cases()*')
                        dose_included_event.remove_property('dose')
                    # minor cleanup: there may still be loose strucs that are assigned Schedule_1 or one of its events as
                    # their Sems (e.g. Then_Chrono). So reassign these to None.
                    for struc in instruction.parse.strucs:
                        if struc.accounted_for_by_sem == schedule_1 or struc.accounted_for_by_sem in schedule_1.events:
                            struc.accounted_for_by_sem = None
                            struc.rules_used.append('*removed_from_sem_by_process_special_duration_cases()*')
                    schedule_1_events = schedule_1.events[:]
                    for event in schedule_1_events:
                        schedule_1.remove_event(event)
                    return versions
        return []

    def combine_periodicity_with_as_needed(version, instruction):
        """ Combine instatoms "Take this medicine every day. Take this medicne as needed" into "Take daily as needed". Same for "use", etc.

        Justification:
        "take 1 tablet. Take this medicine every day. Take this medicine as needed", even though it reflects the meaning of "take 1 tablet daily as needed"
        still has a discombobulating feeling of a contradiction because it isolates "every day" and puts it next to "as needed". The reason the source doesn't
        sound self-contradictory is because "daily" is combined with "as needed" and not juxtaposed with it.

        Unfortunately, we can only combine Periodicity and As_needed at the instruction assembly stage because Periodicity is part of Schedule while all the
        proper properties of Instruction are dealt with separately in process_instruction_stucs_proper(). That happened because we assumed that
        properties such as Indication, As_Needed, See_Instruction should only occur once per instruction and tend to be combined with each other (e.g. as_needed
        with indications). And it is ok if Periodicity is separate from As_Needed if As_needed is combined with an indication (e.g. "take this medicine as needed for headache.")

        So we only catch here "take/use/apply (this medicine)? every day" immediately followed by "take/use this medicine as needed"
        """

        if instruction.as_needed and not instruction.indication and instruction.schedules and instruction.schedules[-1].periodicity:
            pass
        else:
            return

        num_instatoms = len(version.instatoms)
        as_needed = None
        new_instatom = None
        for instatom_index in xrange(-1, -1 * num_instatoms, -1):
            instatom = version.instatoms[instatom_index]
            parse_string = instatom.atom.parses[0].flatten()
            if parse_string in ('DIRECTIVE AS_NEEDED', 'DIRECTIVE ANAPHORA AS_NEEDED'):
                as_needed = instatom.matching_strucs[-1]
                # omit directive - we will take the directive from the periodicity instatom
                break

        if as_needed:
            previous_instatom = version.instatoms[num_instatoms + instatom_index - 1]
            parse_string = previous_instatom.atom.parses[0].flatten()
            if parse_string in ('DIRECTIVE PERIODICITY', 'DIRECTIVE ANAPHORA PERIODICITY'):
                # create strucs_to_match to be the list [directive, periodicity, anaphora, as_needed'] and then try to match it to the dictionary
                directive = previous_instatom.matching_strucs[0]
                periodicity = previous_instatom.matching_strucs[-1]
                strucs_to_match = [directive, periodicity, as_needed]
                candidate_atoms = get_structurally_similar_list_of_atoms_in_dict(strucs_to_match)
                empty_version = Version(version.parse)
                new_versions = try_adding_struc_list_to_versions_list(strucs_to_match, [empty_version], restrict_to_these_atoms_list = candidate_atoms)
                Version.prune_and_sort_versions(new_versions)
                if new_versions and len(new_versions[0].instatoms) == 1 and set([directive, periodicity, as_needed]) <= new_versions[0].matched_strucs:
                    new_instatom = new_versions[0].instatoms[0]

        if new_instatom:
            version.instatoms[num_instatoms + instatom_index - 1:num_instatoms + instatom_index + 1] = [new_instatom]

    def process_instruction(instruction, old_versions):
        """ Extracts the Strucs for the first sentence (core) of the output.

         We want to run try_adding_struc_list_to_versions_list() always with full list of all existion versions, so that we could do checking
        context if we need to. But we only want to return the newly created versions, because we sometimes only want to check
        if the version is addable without committing to adding the new strucs to it (e.g. when we watn to add "then" only if the schedule has substance.).

        """

        parse = instruction.parse

        new_versions = []  # These are versions that don't incorporate old_versions content. We don't add everything to old_versions in order to be able to create sub-versions out of order.
        special_versions = process_special_duration_cases(instruction)
        new_versions += special_versions

        for sched_num, schedule in enumerate(instruction.schedules):
            if sched_num == 0:
                new_schedule_versions = process_schedule(schedule, is_initial_schedule = True)
                if not new_schedule_versions:
                    # If we can't process the first schedule, we should not start sticking "Then blah" after nothing. Just skip to processing Indications (and other stuff at the Instruction level)
                    break
            else:
                if schedule.taper and not (schedule.taper.then_flag and (instruction.schedules[0].duration or instruction.schedules[0].calendar_event)):
                    # For Tapering cases, omit "AFTER THAT:" if there is no duration in the previous schedule. E.g.
                    # "take one capsule at bedtime increase by 1 capsule every day as needed as directed for pain" should not have "AFTER That:"
                    is_initial_schedule = True
                else:
                    is_initial_schedule = False
                new_schedule_versions = process_schedule(schedule, is_initial_schedule = is_initial_schedule)
                # when there are more than 1 schedules, the number of versions may grow exponentially so set max.
                if len(instruction.schedules) == 2:
                    max_number_versions = 4
                else:
                    max_number_versions = 2
                new_schedule_versions.sort(key = lambda version: version.quality_score(), reverse = True)
                new_schedule_versions = new_schedule_versions[:max_number_versions]

            for version in new_schedule_versions:
                version.reorder_instatoms_in_schedule()

            new_versions = Version.concatenate_versions(new_versions, new_schedule_versions)

        instruction_versions_proper = process_instruction_stucs_proper(instruction)
        new_versions = Version.concatenate_versions(new_versions, instruction_versions_proper)
        if instruction_versions_proper:
            for version in new_versions:
                version.reorder_instatoms_in_instruction()
                combine_periodicity_with_as_needed(version, instruction)
        return new_versions



    all_versions = []

    for parse in sig.parses:
        if not sig.raw_sig:
            break
        if not parse.instructions:
            raise UnparseableSig('Parse has no instructions', sig.raw_sig)

        if not parse.instructions[0].is_drugadmin_instruction():
            raise UnparseableSig('First Instruction  of the Sig is not of type DrugAdmin', sig.raw_sig)

        potentially_applicable_atoms_for_whole_sig = prune_out_incompatible_atoms_by_keys(parse)

        initial_version = Version(parse)
        extendable_versions = [initial_version]

        # Try to see if we have a match to the sig as it is:
        tentative_struc_list = parse.strucs
        new_versions = try_adding_struc_list_to_versions_list(tentative_struc_list, extendable_versions)

        # Because these are maximal versions, we don't add them to extendable_versions but only to versions_this_parse.
        versions_this_parse = new_versions

        for instruction in parse.instructions:
            new_versions = process_instruction(instruction, old_versions = extendable_versions)
            extendable_versions = new_versions
            versions_this_parse += extendable_versions
        all_versions += versions_this_parse
        #finished cycling through this parse

    unique_versions = Version.prune_and_sort_versions(all_versions)
    for version in unique_versions:
        version.cleanup_remove_debris()
    sig.versions = unique_versions


    ######## End break_sig_into_atoms()  ##########




#################   Dictionaries and Pickles    ######


class AtomicDictionary(object):

    def __init__(self):
        self.sigs = []
        self.match_labels_list_2_atoms_list = defaultdict(list)

    def get_atom_sig_list(self):
        return self.sigs

    def get_dictionary_size(self):
        sigs = self.get_atom_sig_list()
        return len(sigs)

    def get_sig_from_english_string(self, english_string):
        for sig in self.get_atom_sig_list():
            if sig.raw_sig == english_string:
                return sig
        return None

    def get_string_dictionary(self, english_string):
        sig = self.get_sig_from_english_string(english_string)
        if sig is None:
            return None
        else:
            return sig.dictionary

    def is_string_in_dict(self, english_string):
        return self.get_string_dictionary(english_string) is not None

    def add_new_english_string(self, english_string):
        if self.is_string_in_dict(english_string):
            raise Exception('This english string is already in the atomic dictionary %s' % english_string)
        sig = parse_sig(english_string, debugging_flag = False)
        sig.dictionary = {en: english_string}
        sig.key_set = get_key_set(sig.parses[0].strucs)
        struc_list = reorder_strucs_in_canonical_order_and_omit_punctuation(sig.parses[0].strucs)
        struc_list_string = get_struc_labels(struc_list, delimiter = '|', omit_spaces_and_punctuation = False)
        self.match_labels_list_2_atoms_list[struc_list_string].append(sig)
        self.sigs.append(sig)

    def update_translation_for_locale(self, english_string, locale, translation):
        sig = self.get_sig_from_english_string(english_string)
        if sig is None:
            raise Exception('This english_string -->%s<-- is not in the dictionary at all' % english_string)
        else:
            sig.dictionary[locale] = translation

    def clean_whitespace_from_atoms_dictionary(self):
        """Strip extraneous whitespace (carriage return "\r") frome keys and values in the atoms_dictionary

        Use this function once after loading atoms_dictionary from pickle.
        There is a chance that this would be unnecessary if the pickle file were opened in binary mode using pickle (not cPickle) """

        sigs = self.get_atom_sig_list()
        for sig in sigs:
            new_dict = dict()
            old_dict = sig.dictionary
            for locale, translation in old_dict.items():
                locale = locale.strip()
                translation = translation.strip()
                new_dict[locale] = translation
            sig.dictionary = new_dict


atoms_dictionary = AtomicDictionary()
atoms_dictionary_path = 'transduction/pickles/atoms_dictionary.pkl'
tentative_atomic_dictionary = AtomicDictionary()
tentative_atomic_dictionary_path =  'pickles/tentative_atomic_dictionary.pkl'
approved_transductions_by_raw_sig = {}
approved_transductions_by_raw_sig_path = 'pickles/approved_transductions.pkl'


def unpickle_object(input_file_name):
    input_file = open(input_file_name, 'r')
    object_name = cPickle.load(input_file)
    input_file.close()
    return object_name

def unpickle_atoms_dictionary(pickle_path = None):
    if not pickle_path:
        pickle_path = atoms_dictionary_path
    try:
        dictionary = unpickle_object(pickle_path)
    except IOError as e:
        print('Exception (tolerable only if we are creating a pickled dictionary for the first time): cant open pickle at |%s|. Error: %s' % (pickle_path, e))
        dictionary = AtomicDictionary()
        return dictionary
    dictionary.clean_whitespace_from_atoms_dictionary()
    return dictionary

def get_atoms_dictionary():
    global atoms_dictionary, tentative_atomic_dictionary

    if atoms_dictionary.get_dictionary_size() == 0:
        if debug:
            print('******* Unpickling atoms_dictionary from path %s' % atoms_dictionary_path)
        atoms_dictionary = unpickle_atoms_dictionary(atoms_dictionary_path)
        if debug:
            print('******* Unpickled atoms_dictionary of size %d\n*********\n' % (atoms_dictionary.get_dictionary_size()))

    # When testing, also open up tentative_atomic_dictionary and merge it into atoms_dictionary
    if debug and tentative_atomic_dictionary.get_dictionary_size() == 0:
        print('******* Unpickling tentative_atoms_dictionary from path %s' % tentative_atomic_dictionary_path)
        tentative_atomic_dictionary = unpickle_atoms_dictionary(tentative_atomic_dictionary_path)
        print('******* Unpickled tentative_atomic_dictionary of size %d\n*********\n' % (tentative_atomic_dictionary.get_dictionary_size()))
        print('\n Size of atoms_dict before merger: %d with %d frames' % (atoms_dictionary.get_dictionary_size(), len(atoms_dictionary.match_labels_list_2_atoms_list)))
        atoms_dictionary.sigs+= tentative_atomic_dictionary.sigs
        for (struc_list_string, sig_list) in tentative_atomic_dictionary.match_labels_list_2_atoms_list.items():
            atoms_dictionary.match_labels_list_2_atoms_list[struc_list_string].extend(sig_list)
        print('\n After merging atoms_dict with tent atoms dict, size of atoms_dict %d with %d frames' % (atoms_dictionary.get_dictionary_size(), len(atoms_dictionary.match_labels_list_2_atoms_list)))

def show_atom(raw_sig):
    """ For debugging only.
    """

    sig = atoms_dictionary.get_sig_from_english_string(raw_sig)
    source = 'atoms'
    if not sig:
        sig = tentative_atomic_dictionary.get_sig_from_english_string(raw_sig)
        source = 'tentative_atoms'
    if not sig:
        print('Not found in either atoms_dictionary or tentative_atomic_dictionary')
        return
    print('Raw sig: ->%s<-\nKey Set: %s' % (sig.raw_sig, ', '.join([key for key in sorted(list(sig.key_set))])))
    print sig.show()
    return sig


#################   End Dictionaries and Pickles #####



def test(raw_sigs = None, debugging_flag = None, details_to_show = None, limit_to_incompletely_parsed = None, return_sigs = None, print_developer_msgs = None, do_not_execute = None):
    """
    Runs transduce on raw_sigs corpora. By default, it's Target 1K.
    details_to_show                 Integer. If 2 shows full sig with parse results, if 0, shows only raw_sig
    limit_to_incompletely_parsed    True/False. If True, shows only sigs not fully transduced.
    return_sigs                     True/False. If True and limit_to_incompletely_parsed=True, returns incompletely parsed sig objects.
                                                If True and not limit_to_incompletely_parsed, retursn all sig objects.
    """

    def readfile(fn = None):
        if not fn:
            fn = '../data/all_1172602_unique_occurences_from_Target_Supervalue_and_Logs_thru_1_27_2012.txt'

        print('About to read file %s' % fn)
        fi = open(fn, 'r')
        lines = fi.readlines()
        fi.close()
        print('Read file %s with %d lines' % (fn, len(lines)))
        lines = [trim(line) for line in lines]
        return lines

    if do_not_execute:
        return None

    if raw_sigs == None:
        fn = '../data/959_SIGs.txt'
        raw_sigs = readfile(fn)
    elif type(raw_sigs) == str and raw_sigs[-4:] == '.txt':
        raw_sigs = readfile(raw_sigs)

    start_time = datetime.datetime.now()
    sigs_to_show = []
    incompletely_parsed = 0
    incompletely_segmented = 0
    for i, raw_sig in enumerate(raw_sigs):
        sig = transduce(raw_sig, debugging_flag = debugging_flag)
        is_fully_parsed = False
        is_fully_sem_segmented = False
        for parse in sig.parses:
            if parse.is_fully_parsed():
                is_fully_parsed = True
                if parse.is_fully_sem_segmented():
                    is_fully_sem_segmented = True
                break

        if not is_fully_parsed:
            incompletely_parsed += 1
        if not is_fully_sem_segmented:
            incompletely_segmented += 1

        if not limit_to_incompletely_parsed:
            sigs_to_show.append((sig, i))
        elif len(raw_sigs) == 1 or not is_fully_parsed or not is_fully_sem_segmented:
            sigs_to_show.append((sig, i))

    end_time = datetime.datetime.now()
    time_took = end_time - start_time

    print('Out of %d raw sigs, not fully parsed are %d and parsed but not fully segmented are %d' % (len(raw_sigs), incompletely_parsed, incompletely_segmented))


    for i, (sig, sig_num) in enumerate(sigs_to_show):
        if details_to_show:
            print('\n-------------Sig No. %d -----------------------------------------------' % sig_num)
            sig_descr = sig.pprint(details_to_show = details_to_show)
            print sig_descr
        else:
            print('%4d  %4d  %s' % (i, sig_num, sig.raw_sig))

    print('Out of %d raw sigs, not fully parsed are %d' % (len(raw_sigs), incompletely_parsed))
    #print('transduce for %d sigs took: %d min %02d secs, or %.3f seconds per sig' % (len(raw_sigs), time_took.seconds // 60, time_took.seconds % 60, time_took.total_seconds()/len(raw_sigs)))

    if print_developer_msgs:
        DeveloperMessages.print_developer_messages()

    if return_sigs:
        sigs = [sig for (sig, x) in sigs_to_show]
        return sigs




###### Main

def parse_sig(raw_sig, debugging_flag = None):

    global debug

    # debug is a module variable that if True raises exceptions in case of constraint violations and helps build
    # auxiliary information about the training corpus.

    if debugging_flag: debug = True
    else: debug = False

    get_atoms_dictionary()

    sig = cleanup_raw_sig(raw_sig)
    apply_preprocess_rules(sig)
    apply_struc_identification_rules(sig)
    create_sems(sig)
    modify_sems(sig)
    return sig

def transduce(raw_sig, debugging_flag = None):
    """    Transduces the sig. Returns a list of English transductions (paraphrases) in the decreasing order of likelyhood.    """

    reset_list(errors_list)

    sig = parse_sig(raw_sig, debugging_flag)

    try:
        break_sig_into_atoms(sig)
    except UnparseableSig as e:
        if debugging_flag:
            error_msg = ('\n-------->>>UnparseableSig Error: %s' % e.__str__())
            msg = DeveloperMessages(error_msg, msg_type = 'Error', error_level = 1)
            errors_list.append(msg)
            raise UnparseableSig(error_msg, raw_sig)




    if debug and errors_list:
        msg = DeveloperMessages(raw_sig, msg_type = 'raw_sig', error_level = 0)
        errors_list.insert(0, msg)
        collective_errors_all_sigs.append(errors_list)


    return sig


def process_one_sig(raw_sig, max_number_of_versions = None):
    """ Transduce the raw_sig and return a dictionary mapping locales to lists of translated versions.

        For example, the dictionary may map 'ru_RU' locale to a list of two strings, each string representing a different version of a translation
        into Russian of the same raw_sig. The different translations actually correspond to different versions of the transduction attempt
        of the raw_sig into the fixed set of dictionary atoms.

    max_number_of_versions      Maximum number of versions of a translation/transduction that should be provided.
    """

    if not max_number_of_versions:
        max_number_of_versions = 4

    sig = transduce(raw_sig, debugging_flag = False)
    locale_2_list_of_versions = defaultdict(list)
    for version in sig.versions[:max_number_of_versions]:
        locale_2_translation = version.get_instantiated_translations()
        for (locale, translation) in locale_2_translation.items():
            locale_2_list_of_versions[locale].append(translation)


    return locale_2_list_of_versions

def transduce_sig(raw_sig):
    """
    Return a dict whose keys are locales and whose values are ordered lists of translations of the English paraphrase of raw_sig.

    Return an empty dict if no English paraphrase can be transduced for raw_sig; i.e., if the result is not empty,
    it must include a correct value for the locale 'en_US'. The locale 'en_US' maps to the English paraphrase of raw_sig
    generated by transduction. Other locales map to the respective translation of the English paraphrase.

    Locales are identified by the standard codes, e.g. 'zh_CN' for "Chinese (Simplified)", 'fr_FR' for "French".
    Recognized locale codes, maps between language name and locale code and vice versa, and methods for adding
    new language codes are available through the Locale_Map class, which stores/retrieves them through a pickled dict.

    Translation values are unicode objects which, for printing purposes, should be encoded to utf8 via utf_encode(string) function.
    """

    try:
        translations_dictionary = process_one_sig(raw_sig)
        return translations_dictionary
    except: return {}



if __name__ == "__main__":

    raw_sigs_difficult = ['take one tablet daily 4 times per week',
                          'TAKE ONE TABLET 0 TABLETS ON TUESDAYS AND FRIDAY AND 1 TABLET ON OTHER DAYS',
                          'take one tablet by mouth one time daily skipping sundays',
                          'take one tablet by mouth daily continuously as directed by prescriber. skip placebo.']


    raw_sigs = ['1 qday then bid then tid for two weeks, then 2 qday for three weeks, then 3 qday for four weeks',
                'TAKE TWO TABLETS BY MOUTH ON DAY 1, AND THEN TAKE ONE TABLET BY MOUTH ONCE A DAY FOR DAYS 2-5',
                'Take two tablets on Monday in the morning. Take one tablet on Tuesday in the evening',
                'dissolve 1 tablet in mouth every night at bedtime',
                'take two tablets on monday and wednesday',
                '2 teaspoonfools orally stat then 1 teaspoonful orally 3 times daily until gone',
                'dissolve 1 tablet in mouth every night at bedtime',
                'take one tablet by mouth daily at bedtime for cholesterol',
                'TAKE <<NUM_0>> TABLETS BEFORE BEDTIME.',
                'take 2 tabs now, then for the next 2 weeks take 1 tab every morning and 2 tabs before bedtime for severe pain',
                'take 2 tabs now, then for the next 2 weeks take 1 tab every morning and 2 tabs before bedtime',
                'TAKE 1 TABLET BY MOUTH ONCE DAILY TO TWICE DAILY AS NEEDED FOR SEVERE PAIN',
                'Take one 1 cap daily',
                'Take 1 capsule twice daily with food and water.',
                'TAKE 1 TABLET BY MOUTH ONCE DAILY TO TWICE DAILY AS NEEDED FOR SEVERE PAIN',
                'take 1 tablet by mouth 1 time daily at bedtime',
                'INHALE ONE PUFF BY MOUTH EVERY FOUR HOURS AS NEEDED',
                '1 vial via nebulizer every 4 hours as needed for cough or wheeze',
                'take four tablets by mouth for 2 days decrease by 1 tablet every 2 days then stop',
                'TAKE 2 AND A HALF TABLETS BY MOUTH 2 TIMES A DAY FOR PAIN',
                'swish 1-2 teaspoons four times daily and spit',
                'swish and spite 5ml twice daily',
                'swish 0.5 ounces by mouth twice daily',
                'swish and spite 5ml twice daily',
                'take 2 puffs by mouth twice daily',
                'take one puff by mouth 2to 4 times daily for 30 days as needed for wheeze',
                'take one puff by mouth 1 to 4 times daily for 30 days as needed for wheeze',
                'take 2 puffs by mouth twice daily',
                'TAKE DAILY AS NEEDED.',
                'swish and spit four times daily',
                '2 sprays each nostril for 7 days then 1 spray each nostril daily for 7 days then stop',
                '2 sprays each nostril for 7 days then 1 spray each nostril daily for 7 days then stop',
                '2 sprays each nostril for 7 days then 1 spray each nostril daily for 7 days then stop',
                 'take one capsule sunday and wednesday for infection',         ###
                 'take one capsule by mouth every sunday for 3 months',
                 'inject 10 units at breakfast 12 units at lunch and 10 units at dinner', ######
                 'dissolve 1 tablet sublingually every 5 minutes up to 3 doses as needed then call md',
                 'dissolve 1 tablet sublingually every 5 minutes up to 3 doses as needed then call doctor',
                 'TAPER AS DIRECTED',
                 'taper as directed with food',
                 'INSERT 1 RING VAGINALLY. REMOVE AFTER 3 WEEKS & ALLOW 1WEEK REST.',
                 'take one capsule by mouth every morning take one tablet at noon and take two tablets at bedtime',
                 'DISSOLVE 1 TABLET UNDER THE TONGUE FOR CHEST PAIN MAY REPEAT EVERY 5 MINUTES UP TO 3 TIMES IF NOT RESOLVED CALL 911',
                 'take 1.5-2 capsules by mouth every morning for 7 days as directed by your md and 3-4 by mouth every evening',
                 'insert ring in vagina and leave in for 3 weeks then remove for 1 week, then repeat',
                 'place ring in vagina and leave in for 3 weeks then removeWISH AND SWALLOW 30 ML.|USE THIS MEDICINE EVERY 6 HOURS.|USE THIS MEDICINE FOR 1 WEEK(S). for 1 week, then repeat',
                 'APPLY ONE PATCH ONCE WEEKLY FOR 3 WEEKS THEN OFF FOR ONE WEEK']

    fn = '../data/Temp/queries_top_1000_thru_8_30_2012.txt'
    fn = '../data/Temp/dictionary.txt'

    if False:
        start_time = datetime.datetime.now()
        sigs = test(raw_sigs = raw_sigs[0:1], debugging_flag = False, details_to_show = 10, limit_to_incompletely_parsed = False, return_sigs = True, print_developer_msgs = True, do_not_execute = False)
        sig = sigs[0]


        if False:
            sig = transduce(raw_sigs[0], debugging_flag = True)

        parse = sig.parses[0]
        strucs = parse.strucs
        time_took = datetime.datetime.now() - start_time

        print sig.pprint(details_to_show = 0)
        #DeveloperMessages.print_developer_messages()

        #print('transduction took %.3f seconds' % (time_took.total_seconds()))
        if False:
            translations_dictionary = transduce_sig(raw_sigs[0])
            for (locale, translations) in translations_dictionary.items():
                print locale
                for translation in translations:
                    print('  %s' % translation)

    if False:
        show_atom('USE THIS MEDICINE FOR <<NUM_0>>-<<NUM_1>> WEEKS.')


    if False:
        txt1 = 'Take 1 tablet by mouth every day'
        txt2 = 'take 1.5 tabs by mouth 2 to 4 times daily for 30 days as needed for wheezing'
        txt3 = 'take 5 tabs by mouth 3 to 4 times daily for 30 days as needed for'
        translations_dictionary = transduce_sig(txt2)
        #print translations_dictionary
        #print 'Just printed result of transduce_sig(txt2)'
        for (locale, translations) in translations_dictionary.items():
            print locale
            for translation in translations:
                print('  %s' % translation)




