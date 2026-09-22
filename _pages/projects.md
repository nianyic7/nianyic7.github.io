---
title: "Research"
permalink: /projects/
classes: wide
---

{% assign themes = site.data.research %}

<nav class="research-themes" aria-label="Research themes">
  <ul>
  {% for theme in themes %}
    <li><a href="#{{ theme.id }}">{{ theme.title }}</a></li>
  {% endfor %}
  </ul>
</nav>

{% for theme in themes %}
<section class="research-theme" id="{{ theme.id }}">
  <h2>{{ theme.title }}</h2>
  <div class="research-theme__intro{% if theme.image %} has-image{% endif %}">
    <p>{{ theme.summary }}</p>
    {% if theme.image %}<img src="{{ theme.image | relative_url }}" alt="{{ theme.title }}">{% endif %}
  </div>
  <h3>Selected papers</h3>
  <ol class="research-theme__papers">
  {% for p in theme.papers %}
    <li>
      <span class="paper-title">{{ p.title }}</span><br>
      <span class="paper-meta">{% if p.me %}<strong>{{ p.authors }}</strong>{% else %}{{ p.authors }}{% endif %} <em>{{ p.venue }}</em> ({{ p.year }})
      · <a href="https://doi.org/{{ p.doi }}">DOI</a>{% if p.arxiv %} · <a href="https://arxiv.org/abs/{{ p.arxiv }}">arXiv</a>{% endif %}</span>
    </li>
  {% endfor %}
  </ol>
</section>
{% endfor %}
