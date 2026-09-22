---
title: "Miscellaneous"
permalink: /misc/
---

## Notes for collaborators and myself

{% assign md_notes = site.notes | sort: "date" | reverse %}
{% assign pdf_notes = site.static_files | where_exp: "f", "f.path contains '/assets/notes/'" | where: "extname", ".pdf" %}

{% if md_notes.size == 0 and pdf_notes.size == 0 %}
*Nothing here yet.*
{% else %}
<ul class="notes-list">
{% for note in md_notes %}
  <li>
    <a href="{{ note.url | relative_url }}">{{ note.title | default: note.slug }}</a>
    {% if note.date %}<span class="notes-list__meta">· {{ note.date | date: "%b %-d, %Y" }}</span>{% endif %}
    {% if note.description %}<br><span class="notes-list__desc">{{ note.description }}</span>{% endif %}
  </li>
{% endfor %}
{% for f in pdf_notes %}
  <li>
    <a href="{{ f.path | relative_url }}">{{ f.basename | replace: "_", " " | replace: "-", " " }}</a>
    <span class="notes-list__meta">· PDF</span>
  </li>
{% endfor %}
</ul>
{% endif %}
