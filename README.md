# nianyic7.github.io

Personal website of Nianyi Chen, built with Jekyll and the
[Minimal Mistakes](https://github.com/mmistakes/minimal-mistakes) theme (4.28.1).

## Local preview

```sh
bundle install
bundle exec jekyll serve
```

Then open http://localhost:4000.

## Where things live

- `_pages/about.md`: home page
- `_pages/publications.md`: publications page (links to Google Scholar)
- `_pages/cv.md` + `_data/cv.yml` + `assets/pdf/CV_Nianyi.pdf`: CV (hidden for now: `published: false`, and commented out in the nav and sidebar)
- `_pages/misc.md`: Miscellaneous page; notes are listed automatically from `_notes/` (Markdown) and `assets/notes/` (PDF)
- `_pages/projects.md` + `_data/research.yml`: Research page (themes and selected papers)
- `_data/navigation.yml`: top navigation
- `_config.yml`: site settings and sidebar author profile

Pushing to `master` builds the site with GitHub Actions and deploys it to the `gh-pages` branch.

## Adding notes for collaborators

Notes appear automatically under "Notes for collaborators and myself" on the Miscellaneous page.

- **Markdown:** add `_notes/my-note.md` starting with front matter such as
  ```yaml
  ---
  title: "My note"
  date: 2026-09-22
  description: "Optional one-line summary"
  ---
  ```
  It is published at `/notes/my-note/`. LaTeX math (`$...$`) works.
- **PDF:** drop `my_note.pdf` into `assets/notes/`. It is published at `/assets/notes/my_note.pdf`.

Everything here is public: anyone with the link (or who browses the Miscellaneous page) can open it.
