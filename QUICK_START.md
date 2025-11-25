# Quick Start Guide

Get your thesis portfolio live in 15 minutes.

## What You Just Got

A complete GitHub Pages website with:
- ✓ Dark academic theme
- ✓ Chapter-based organization matching your thesis structure
- ✓ Templates for Chapters 3 & 4 (Campylobacter and Salmonella)
- ✓ Reusable methodology templates
- ✓ Example scripts with documentation
- ✓ CSS and HTML already configured
- ✓ Private by default (public when you're ready)

## The 5-Minute Setup

### 1. Create a GitHub Repository

```bash
# On GitHub.com:
# 1. Click + → New repository
# 2. Name: thesis-portfolio
# 3. Visibility: Private
# 4. Create repository
```

### 2. Push Your Local Copy

```bash
cd thesis-portfolio
git init
git add .
git commit -m "Initial commit: thesis portfolio"
git remote add origin https://github.com/YOUR_USERNAME/thesis-portfolio.git
git branch -M main
git push -u origin main
```

### 3. Enable GitHub Pages

```bash
# On GitHub.com:
# 1. Settings → Pages
# 2. Source: main branch, / (root)
# 3. Save
# Wait 1-2 minutes for build to complete
```

### 4. View Your Site

Visit: `https://YOUR_USERNAME.github.io/thesis-portfolio/`

### 5. Share with Reviewers

```bash
# On GitHub.com:
# 1. Settings → Collaborators
# 2. Add people with GitHub usernames
# 3. Share link from step 4
```

Done! Your site is live.

---

## Now What?

### Quick Wins (30 minutes)

1. **Edit `about.md`**
   - Add your institution name
   - Update research focus
   - Add contact info

2. **Update `_config.yml`**
   - Change `title` to your name
   - Update `author` and `email`
   - Optionally add a logo

3. **Customize colors** (optional)
   - Edit `assets/css/style.css`
   - Change `--accent: #4a9eff;` to your favorite color

### Add Your Content (1-2 hours)

1. **For Chapter 3 (Campylobacter):**
   ```
   chapters/03-campylobacter-plasmids/
   ├── 3.2-methodology/
   │   └── 3.2.1-assembly/
   │       ├── README.md              ← Edit this
   │       ├── assembly_qc.sh         ← Replace with your script
   │       └── assembly_stats.R       ← Replace with your script
   └── 3.4-data/
       └── [Add your data files here]
   ```

2. **For Chapter 4 (Salmonella):**
   - Same structure, different section
   - Follow the templates

3. **Update index.md pages:**
   - `chapters/03-campylobacter-plasmids/index.md`
   - `chapters/04-salmonella-surveillance/index.md`
   - Add real file paths, descriptions, and links

### Workflow for Adding Files

```bash
# Add a script
cp my_analysis.sh chapters/03-campylobacter-plasmids/3.2-methodology/3.2.1-assembly/

# Add data
cp my_data.csv chapters/03-campylobacter-plasmids/3.4-data/processed/

# Update the chapter index.md to list your files

# Commit and push
git add chapters/03-campylobacter-plasmids/
git commit -m "Add assembly scripts and data for chapter 3"
git push
```

Site updates automatically within 1-2 minutes.

---

## File Location Guide

**Where your Campylobacter files go:**
```
chapters/03-campylobacter-plasmids/
├── 3.2-methodology/          Scripts and workflows
├── 3.3-results/              Output files, tables, figures
└── 3.4-data/                 Raw and processed data
```

**Where your Salmonella files go:**
```
chapters/04-salmonella-surveillance/
├── 4.2-methodology/          Scripts and workflows
├── 4.3-results/              Output files, tables, figures
└── 4.4-data/                 Raw and processed data
```

**How reviewers navigate:**
1. Home page → Overview
2. Click Chapters → See all chapters
3. Click Chapter 3 → See all subsections
4. View methodology, results, data
5. Click on script → View code in browser or download
6. Click on Excel file → Download and open

---

## Common Tasks

### View Your Site Locally

```bash
# Install Jekyll (one time)
gem install bundler jekyll

# Build and serve
bundle exec jekyll serve

# Visit http://localhost:4000/
```

### Add a New Chapter

```bash
# 1. Create folder
mkdir -p chapters/05-chapter-name

# 2. Create index.md
cp chapters/03-campylobacter-plasmids/index.md chapters/05-chapter-name/

# 3. Edit the new index.md with your content

# 4. Add link in chapters.md

# 5. Commit and push
git add chapters/05-chapter-name/
git commit -m "Add Chapter 5: [Title]"
git push
```

### Update a Script

```bash
# 1. Replace the file
cp my_updated_script.sh chapters/03-campylobacter-plasmids/3.2-methodology/3.2.1-assembly/

# 2. Update the README if methods changed

# 3. Commit and push
git add chapters/03-campylobacter-plasmids/3.2-methodology/3.2.1-assembly/
git commit -m "Update assembly QC script"
git push
```

### Make Site Public (When Ready)

```bash
# On GitHub.com:
# 1. Settings → General
# 2. Danger zone → Change visibility → Public
# 3. Confirm

# Your site stays at same URL, now everyone can see it
```

---

## Templates Available

We've included templates to make things easier:

- **`METHODOLOGY_TEMPLATE.md`** - Use for documenting any method section
- **`example_analysis.sh`** - Use as template for bash scripts
- **Chapter templates** - Chapters 3 & 4 show the structure to follow

Just copy and modify these for your needs.

---

## Helpful Resources

**For GitHub Pages:**
- [GitHub Pages docs](https://docs.github.com/en/pages)
- [Jekyll markdown guide](https://kramdown.gettalong.org/quickref.html)

**For your workflow:**
- Edit files in your code editor
- Test locally with `jekyll serve`
- Push when ready: `git push`
- Site updates automatically

**For styling:**
- Main CSS: `assets/css/style.css`
- Site config: `_config.yml`
- Page template: `_layouts/default.html`

---

## Troubleshooting

**"Page not found" when I visit the site?**
- Site is building (takes 1-2 min). Check again.
- Go to Actions tab to see build status.

**Content not showing?**
- Check you added front matter to .md files:
  ```yaml
  ---
  layout: default
  title: Page Title
  ---
  ```
- Hard refresh: Ctrl+Shift+R

**Can't push to GitHub?**
- Run `git pull origin main` first
- Then `git push origin main`

**Local jekyll serve not working?**
- Make sure Ruby is installed: `ruby --version`
- Install gems: `bundle install`
- Run: `bundle exec jekyll serve`

---

## Next Steps

1. ✅ Follow the 5-minute setup above
2. ✅ View your site at `https://YOUR_USERNAME.github.io/thesis-portfolio/`
3. ✅ Edit `about.md` and `_config.yml`
4. → Add your actual scripts and data
5. → Update chapter index.md files
6. → Share with reviewers
7. → Make public when ready for defense

---

## Questions?

- **Jekyll issues?** Check [jekyllrb.com](https://jekyllrb.com)
- **GitHub issues?** Check [GitHub help](https://docs.github.com)
- **Your setup?** Review [GITHUB_SETUP.md](GITHUB_SETUP.md) for detailed instructions

---

**You've got this! Your reviewers are going to be impressed with how organized everything is.** 🎓

---

*Last updated: November 2025*
