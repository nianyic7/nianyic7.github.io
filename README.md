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
- `_pages/publications.md` + `_bibliography/papers.bib`: publications list (via jekyll-scholar, template in `_layouts/bib.html`)
- `_pages/cv.md` + `_data/cv.yml` + `assets/pdf/CV_Nianyi.pdf`: CV
- `_projects/`: project pages
- `_data/navigation.yml`: top navigation
- `_config.yml`: site settings and sidebar author profile

Pushing to `master` builds the site with GitHub Actions and deploys it to the `gh-pages` branch.
