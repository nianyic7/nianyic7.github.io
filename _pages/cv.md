---
title: "CV"
permalink: /cv/
---

[Download my full CV (PDF)]({{ "/assets/pdf/CV_Nianyi.pdf" | relative_url }}){: .btn .btn--primary}

{% for section in site.data.cv %}
## {{ section.title }}

{% if section.type == "map" %}
{% for item in section.contents %}
- **{{ item.name }}:** {{ item.value }}
{% endfor %}
{% elsif section.type == "time_table" %}
{% for item in section.contents %}
**{{ item.title }}**, {{ item.institution }} ({{ item.year }})
{% for d in item.description %}{% if d.title %}
- {{ d.title }}{% for sub in d.contents %}
  - {{ sub }}{% endfor %}{% else %}
- {{ d }}{% endif %}{% endfor %}

{% endfor %}
{% elsif section.type == "list" %}
{% for item in section.contents %}
- {{ item }}
{% endfor %}
{% else %}
{{ section.contents }}
{% endif %}
{% endfor %}
