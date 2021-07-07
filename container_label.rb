class ContainerLabel < ActiveRecord::Base

  belongs_to :drug
  belongs_to :sig
  belongs_to :paraphrase
  belongs_to :locale
  belongs_to :user
  belongs_to :location
  has_one :patient_education_sheet
  has_one :non_standard_request
  has_many :container_label_warnings, :order => 'position', :dependent => :destroy
  has_one  :corporation, :through => :location

  LEGAL_SERIALIZATION_ATTRS = [ :id, :drug_id, :locale_id, :custom, :reference,
                                :created_at, :sig_token, :paraphrase_token,
                                :translation_token, :translation_token,
                                :translation_image, :container_label_warnings, :label_name ]

  accepts_nested_attributes_for :container_label_warnings

  validates_presence_of :drug_name, :locale, :reference, :pharmacist_name
  validates_presence_of :drug, :message => "verify the drug is not listed", :if => Proc.new{|label| label.drug_not_listed.to_s != "1" }
  validates_acceptance_of :charge_for_custom_translation, :message => "You must indicate your understanding.", :on => :create, :if => Proc.new{|label| label.custom? }
  validate :validates_can_make_custom_translation_requests, :if => Proc.new{|label| label.custom?}
  validates_presence_of :user

  before_validation_on_create :initialize_container_label_warnings, :backfill_drug_name
  before_create :populate_transient_fields
  after_create :build_and_save_pdf

  has_attached_file(*PaperclipAttachmentOptions.new(:label_pdf, self).build)

  attr_accessor :drug_not_listed

  delegate :setting, :to => :corporation_id

  named_scope :created_since, lambda{|created_since|
    { :conditions  => [ "container_labels.created_at >= ?", created_since.to_datetime ] }
  }

  named_scope :today, {
    :conditions => [ "container_labels.created_at >= ?", Date.today ]
  }

  named_scope :with_corporation, lambda { |corporation|
    {
      :joins => [:location],
      :conditions => [ "locations.corporation_id= ?", corporation.id ]
    }
  }

  named_scope :select_location_locale_and_created_at_date, {
    :select => "container_labels.location_id, container_labels.locale_id, container_labels.created_at",
    :conditions => "container_labels.location_id IS NOT NULL AND container_labels.locale_id IS NOT NULL"
  }

  named_scope :order_by, lambda{|order|
    { :order => order }
  }

  named_scope :custom, {
    :conditions => { :custom => true }
  }

  named_scope :standard, {
    :conditions => { :custom => false }
  }

  named_scope :with_user_visibility, lambda{|user|
    if user.is_a? ApiUser
      if user.is? :api_admin
        #:include => [:corporation, {:corporation => :provider}]
        {
          :joins => [:location, {:location => [:corporation, {:corporation => :provider}]}],
          :conditions => [ "container_labels.location_id=locations.id AND locations.corporation_id=corporations.id AND corporations.provider_id=#{user.provider.id}" ]
        }
      else
        { :conditions => { :user_id => user } }
      end
    else
      if user.available_locations.empty?
        { :conditions => "1=2" }
      else
        # TODO: doesn't take into account that corp_admin should be able to see all for corporation
        { :conditions => [ "location_id IN (#{user.available_locations.map(&:id).join(",")})" ] }
      end
    end
  }

  include AASM
  aasm_column :state
  aasm_initial_state :undownloaded
  aasm_state :undownloaded
  aasm_state :downloaded
  aasm_event :download do
    transitions :to => :downloaded, :from => :undownloaded
  end

  def standard?
    !self.custom?
  end

  def container_label_warnings_to_params
    ContainerLabelWarning.to_params(self.container_label_warnings)
  end

  def build_and_save_pdf
    output = StringIO.new
    if self.location.corporation.name.downcase.include?('walgreens')
      pdf_template = IO.read("#{RAILS_ROOT}/app/views/container_labels/walgreen.pdf.pdfbuilder")
    else
      pdf_template = IO.read("#{RAILS_ROOT}/app/views/container_labels/show.pdf.pdfbuilder")
    end

    pdf = PDF::Wrapper.new(output, :paper => :LETTER, :margin_bottom => 0)
    pdf.font locale.font_name
    pdf.instance_variable_set(:@container_label, self)
    pdf.instance_eval do
      eval pdf_template
    end
    pdf.finish
    self.label_pdf = Tempfile.new(pdf_file_name).tap { |f| f.write(output.string) }
    save
  end

  def parse_pdf_template
    html_template = IO.read("#{RAILS_ROOT}/app/views/container_labels/show.pdf.haml")
    haml = Haml::Engine.new(html_template)
    haml.render(self, :label => self)
  end

  def pdf_file_name
    "#{id}-#{reference}.pdf"
  end

  def create_pdf
    build_and_save_pdf unless custom?
  end

  def self.per_page
    Babelscrip::PER_PAGE
  end

  def deliver_notification_non_standard_request_new
    if non_standard_request && non_standard_request.new?
      NotificationMailer.deliver_non_standard_request_new(self.non_standard_request)
    else
      raise NotificationError
    end
  end

  def token
    translation_token
  end

  def to_xml(options={})
    attrs = serialization_attrs_from_options(options)
    image_options = options.only :paper, :font_size, :orientation

    options[:indent] ||= 2
    xml = options[:builder] ||= Builder::XmlMarkup.new(:indent => options[:indent])
    xml.instruct! unless options[:skip_instruct]

    xml.tag! 'container-label' do
      xml.tag! 'id',                id, :type => :integer                    if attrs.include? :id
      xml.tag! 'drug-id',           drug.try(:id), :type => :integer         if attrs.include? :drug_id
      xml.tag! 'locale-id',         locale.id, :type => :integer             if attrs.include? :locale_id
      xml.tag! 'custom',            custom, :type => :boolean                if attrs.include? :custom
      xml.tag! 'reference',         reference                                if attrs.include? :reference
      xml.tag! 'created-at',        created_at.xmlschema, :type => :datetime if attrs.include? :created_at
      xml.tag! 'sig-token',         sig_token                                if attrs.include? :sig_token
      xml.tag! 'paraphrase-token',  paraphrase_token                         if attrs.include? :paraphrase_token
      xml.tag! 'translation-token', translation_token                        if attrs.include? :translation_token

      if attrs.include? :translation_image
        image_format = (options[:image_format] || "png").downcase
        content = LabelMaker::Base64.(self, image_options.stringify_options, image_format).render

        xml.tag! 'translation-image', content, type: "image/#{image_format}"
      end

      if attrs.include? :container_label_warnings
        xml.tag! 'container-label-warnings', :type => :array do
          container_label_warnings.each do |clw|
            xml.tag! 'container-label-warning' do
              xml.tag! 'position',                  clw.position, :type => :integer
              xml.tag! 'warning-token',             clw.warning.token
              xml.tag! 'warning-translation-token', clw.warning_translation.token
              if attrs.include? :translation_image
                image_format = (options[:image_format] || "png").downcase
                content = LabelMaker::Base64.(clw, image_options.stringify_options, image_format).render

                xml.tag! 'translation-image', content, type: "image/#{image_format}"
              end
            end
          end
        end
      end

    end
  rescue
    xml = Builder::XmlMarkup.new
    xml.instruct! unless options[:skip_instruct]
    xml.tag! 'error', $!.message
  end

  def to_json(options={})
    attrs = serialization_attrs_from_options(options)

    cl = ActiveSupport::OrderedHash.new
    cl['id']                = id                if attrs.include? :id
    cl['drug-id']           = drug.try(:id)     if attrs.include? :drug_id
    cl['locale-id']         = locale.id         if attrs.include? :locale_id
    cl['custom']            = custom            if attrs.include? :custom
    cl['reference']         = reference         if attrs.include? :reference
    cl['created-at']        = created_at        if attrs.include? :created_at
    cl['sig-token']         = sig_token         if attrs.include? :sig_token
    cl['paraphrase-token']  = paraphrase_token  if attrs.include? :paraphrase_token
    cl['translation-token'] = translation_token if attrs.include? :translation_token

    if attrs.include? :container_label_warnings
      cl['container-label-warnings'] = container_label_warnings.map do |clw|
        clw_hash = ActiveSupport::OrderedHash.new
        clw_hash['position']                  = clw.position
        clw_hash['warning-token']             = clw.warning.token
        clw_hash['warning-translation-token'] = clw.warning_translation.token
        clw_hash
      end
    end

    {'container-label' => cl}.to_json
  end

  private

  def populate_transient_fields
    # AR Reload
    if self.sig.nil? && self.sig_token
      self.sig = Sig.find_by_token(self.sig_token)
    else
      self.sig.try(:reload)
    end

    if self.paraphrase.nil? && self.paraphrase_id
      self.paraphrase = Paraphrase.find_by_id(self.paraphrase_id)
    else
      self.paraphrase.try(:reload)
    end

    # Copy over
    self.sig_token = self.sig.try(:token)
    self.paraphrase_token = self.paraphrase.try(:token)
    self.translation_token = self.paraphrase.try(:translation_for_locale, self.locale).try(:token)
    # Remove transient relationships
    self.sig, self.paraphrase = nil, nil
  end

  # XXX: Hack until rails 2.3.6
  def initialize_container_label_warnings
    container_label_warnings.each{|clw| clw.container_label = self }
  end

  def serialization_attrs_from_options(options={})
    exclude = options[:exclude]
    exclude = [exclude] unless exclude.is_a?(Array) || exclude.nil?
    only = options[:only]
    only = [only] unless only.is_a?(Array) || only.nil?

    attrs = LEGAL_SERIALIZATION_ATTRS
    attrs -= exclude if exclude
    attrs &= only if only

    attrs
  end

  def validates_can_make_custom_translation_requests
    if custom?
      if not location.corporation.custom_translation_requests_enabled?
        self.errors.add_to_base("Your corporation is not permitted to make custom translation requests.")
      elsif location.above_custom_translation_limit?
        limit = location.non_standard_labels_per_month || 0
        if limit == 0
          self.errors.add_to_base("Your location is not permitted to make custom translation requests.")
        else
          self.errors.add_to_base("Your location is over its quota of #{limit} custom translation requests for per month.")
        end
      end
    end
  end

  def backfill_drug_name
    if drug
      self.drug_name = drug.name
    end
  end

end
