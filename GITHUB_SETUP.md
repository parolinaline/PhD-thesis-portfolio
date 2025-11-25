# GitHub Pages Setup Guide

## Prerequisites

- GitHub account
- Git installed locally
- Basic familiarity with GitHub

## Step 1: Create a GitHub Repository

1. Go to [GitHub](https://github.com) and log in
2. Click the **+** icon in top right → **New repository**
3. Fill in the form:
   - **Repository name**: `thesis-portfolio`
   - **Description**: "PhD Thesis Supplementary Materials"
   - **Visibility**: **Private** (until ready to make public)
   - **Initialize with README**: Leave unchecked (we already have one)
   - **Add .gitignore**: Select "Node" (optional, we have our own)
   - Click **Create repository**

## Step 2: Push Your Local Repository to GitHub

In your local terminal, navigate to the thesis-portfolio directory:

```bash
cd thesis-portfolio

# Initialize git if not already done
git init

# Add all files
git add .

# Create initial commit
git commit -m "Initial commit: thesis portfolio structure with chapters 3 and 4"

# Add the remote repository
git remote add origin https://github.com/YOUR_USERNAME/thesis-portfolio.git

# Rename branch to main (GitHub default)
git branch -M main

# Push to GitHub
git push -u origin main
```

Replace `YOUR_USERNAME` with your actual GitHub username.

## Step 3: Enable GitHub Pages

1. Go to your repository on GitHub: `https://github.com/YOUR_USERNAME/thesis-portfolio`
2. Click **Settings** (right side, top menu)
3. In the left sidebar, click **Pages**
4. Under "Source", select:
   - **Branch**: `main`
   - **Folder**: `/ (root)`
5. Click **Save**
6. GitHub will build your site automatically (takes 1-2 minutes)
7. Once complete, you'll see:
   - "Your site is live at: https://YOUR_USERNAME.github.io/thesis-portfolio/"

## Step 4: Configure Jekyll (Optional but Recommended)

Create a `Gemfile` in your root directory for dependency management:

```ruby
source "https://rubygems.org"

gem "github-pages", "~> 227", group: :jekyll_plugins

group :jekyll_plugins do
  gem "jekyll-feed", "~> 0.17"
  gem "jekyll-seo-tag", "~> 2.8"
  gem "jekyll-remote-theme", "~> 0.4"
end

platforms :mingw, :x64_mingw, :mswin, :jruby do
  gem "tzinfo", ">= 1", "< 3"
  gem "tzinfo-data"
end

gem "wdm", "~> 0.1", :platforms => [:mingw, :x64_mingw, :mswin]
```

Then install and test locally:

```bash
bundle install
bundle exec jekyll serve
# Visit http://localhost:4000/
```

## Step 5: Add Collaborators (for Reviewers)

To share private access with your thesis committee:

1. Go to **Settings** → **Collaborators**
2. Click **Add people**
3. Enter GitHub username or email of reviewer
4. Select permission level (Maintain/View)
5. They'll receive an invitation via email

## Step 6: Customize Your Site (Optional)

### Update Site Title and Description

Edit `_config.yml`:

```yaml
title: "Your Name - PhD Thesis Supplementary Materials"
description: "Comprehensive documentation, scripts, and data for thesis chapters"
author: "Your Full Name"
email: "your.email@institution.ac.nz"
```

### Add Your Logo

Place an image in `assets/images/logo.png` and update `_config.yml`:

```yaml
logo: /assets/images/logo.png
```

### Customize Colors

Edit color scheme in `assets/css/style.css`:

```css
:root {
  --accent: #4a9eff;  /* Change blue to your preference */
  --text-primary: #e0e0e0;
  --bg-primary: #1a1a1a;
}
```

## Step 7: Update About Page

Edit `about.md` with:
- Your actual institution name
- Funding agencies
- Contact information
- Research focus details

## Adding Content

### Add a New Chapter

1. Create folder: `chapters/XX-chapter-name/`
2. Copy structure from `chapters/03-campylobacter-plasmids/` or `chapters/04-salmonella-surveillance/`
3. Create `index.md` with chapter content
4. Add chapter link in `chapters.md`
5. Commit and push:

```bash
git add chapters/XX-chapter-name/
git commit -m "Add Chapter XX: [Title]"
git push
```

### Add Scripts

1. Create subfolder in appropriate chapter/methodology section
2. Add your script file
3. Create a `README.md` explaining the script (use `METHODOLOGY_TEMPLATE.md` as guide)
4. Add links in the chapter's index.md
5. Commit and push:

```bash
git add chapters/XX/XX.Y-section/script.sh chapters/XX/XX.Y-section/README.md
git commit -m "Add analysis script for [subsection]"
git push
```

## Updating the Site

The site rebuilds automatically when you push to GitHub. The process typically takes 1-2 minutes.

To see build status:
1. Go to repository
2. Click **Actions** tab
3. Look for recent "pages build and deployment" workflow
4. Green checkmark = successful

## Local Development Workflow

Recommended workflow for making changes:

```bash
# Create a new branch for your changes
git checkout -b feature/add-chapter-5

# Make your changes (add files, edit markdown, etc.)
nano chapters/05-comparative/index.md
mkdir -p chapters/05-comparative/data

# Test locally
bundle exec jekyll serve
# Visit http://localhost:4000/ to review

# Commit changes
git add chapters/05-comparative/
git commit -m "Add Chapter 5: Comparative Analysis"

# Push to GitHub
git push origin feature/add-chapter-5

# Create Pull Request on GitHub (optional, for review)
# Or simply push to main when ready
git checkout main
git merge feature/add-chapter-5
git push origin main
```

## Troubleshooting

### Site Not Showing

- Check **Actions** tab for build errors
- Review error log details
- Common issues:
  - YAML syntax errors in `_config.yml`
  - Missing required fields in front matter
  - Broken internal links

### Markdown Not Rendering Properly

- Ensure `.md` file has front matter:
```yaml
---
layout: default
title: Page Title
---
```
- Check for unclosed code blocks (triple backticks)
- Validate YAML syntax

### CSS/Styling Not Appearing

- Hard refresh browser: Ctrl+Shift+R (or Cmd+Shift+R on Mac)
- Clear browser cache
- Check CSS file path is correct in HTML

### Git Push Issues

```bash
# If push is rejected:
git pull origin main  # Pull latest changes
git push origin main  # Try again

# If local changes conflict:
git status  # See what changed
git diff    # Review differences
```

## Making Site Public (After Defense)

When ready to make your repository public:

1. Go to **Settings** → **General**
2. Scroll to "Danger zone"
3. Click **Change visibility**
4. Select **Public**
5. Confirm

Your site will remain at the same URL.

## Archiving for Long-term Preservation

Consider these options:

### Zenodo (Free, 50GB limit)
1. Go to [zenodo.org](https://zenodo.org)
2. Connect GitHub account
3. Select your repository
4. Zenodo will create an archive with DOI
5. Get permanent URL and cite in thesis

### Institutional Repository
- Check if your institution has a repository system
- Upload your repository structure there
- Get permanent institutional URL

### GitHub Archive
GitHub automatically preserves all public repositories in the Internet Archive.

## Support

For GitHub Pages issues:
- [GitHub Pages Documentation](https://docs.github.com/en/pages)
- [Jekyll Documentation](https://jekyllrb.com/docs/)
- [GitHub Community Forums](https://github.community/)

For your site specifically:
- Check [Your repo]/Actions for build logs
- Review _config.yml syntax
- Test locally with `jekyll serve`

---

**Next Steps:**
1. ✓ Create GitHub repository
2. ✓ Enable GitHub Pages
3. ✓ Add collaborators for your committee
4. → Start populating your chapters with actual data and scripts
5. → Review site at https://YOUR_USERNAME.github.io/thesis-portfolio/

Good luck with your thesis!
