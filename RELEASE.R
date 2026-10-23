# Steps to release a new version of the EpiEstim package

# Bump version
file.edit("DESCRIPTION")

# Add changes to NEWS
file.edit("NEWS.md")

# Standard checks ------------------------------------------
devtools::test() # Use Ctrl-Shift-T to test with the IDE
devtools::run_examples()
devtools::check() # Ctrl-Shift-E to check with the IDE

# Additional quality checks -------------------------------
lintr::lint_package() # Pre-check lints
spelling::spell_check_package() # Check spelling in docs
urlchecker::url_check() # Check URLs in documentation

# Check docs (as needed)
pkgdown::build_site()

# Platform checks
rhub::rhub_check() # Check on multiple platforms
devtools::check_win_devel() # Check on Windows R-devel

# Push changes
# Once CI clears, request review, merge PR upon approval
# Create GitHub release
