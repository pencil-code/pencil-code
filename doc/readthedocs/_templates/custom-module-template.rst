{{ fullname | escape | underline}}

.. automodule:: {{ fullname }}

   .. rubric:: {{ _('Revision history') }}
   .. git_changelog::
      :fullname: {{ fullname }}

   {% if attributes %}
   .. rubric:: {{ _('List of attributes') }}
   .. autosummary::
   {% for item in attributes %}
      {{ item }}
   {%- endfor %}
   {% endif %}

   {% if exceptions %}
   .. rubric:: {{ _('List of exceptions') }}
   .. autosummary::
   {% for item in exceptions %}
      {{ item }}
   {%- endfor %}
   {% endif %}

   {% if functions %}
   .. rubric:: {{ _('List of functions') }}
   .. autosummary::
   {% for item in functions %}
      {{ item }}
   {%- endfor %}
   {% endif %}

   {% if classes %}
   .. rubric:: {{ _('List of classes') }}
   .. autosummary::
   {% for item in classes %}
      {{ item }}
   {%- endfor %}
   {% endif %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: {{ _('Attributes') }}

   {% for item in attributes %}
   .. autodata:: {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block exceptions %}
   {% if exceptions %}
   .. rubric:: {{ _('Exceptions') }}

   {% for item in exceptions %}
   .. autoexception:: {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block functions %}
   {% if functions %}
   .. rubric:: {{ _('Functions') }}

   {% for item in functions %}
   .. autofunction:: {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block classes %}
   {% if classes %}
   .. rubric:: {{ _('Classes') }}

   {% for item in classes %}
   .. autoclass:: {{ item }}
      :special-members: __init__
      :show-inheritance:
      :members:
      :undoc-members:
   {%- endfor %}
   {% endif %}
   {% endblock %}

{% block modules %}
{% if modules %}
.. rubric:: Modules

.. autosummary::
   :toctree:
   :template: custom-module-template.rst
   :recursive:
{% for item in modules %}
   {{ item }}
{%- endfor %}
{% endif %}
{% endblock %}