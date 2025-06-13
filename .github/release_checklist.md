Release checklist
- [ ] Check outstanding issues on JIRA and GitHub.
- [ ] Check README looks fine.
- [ ] Create a release branch.
  - [ ] Change current development version in `CHANGELOG.rst` to stable version.
- [ ] Merge the release branch into `main`.
- [ ] Created an annotated tag with the stable version number. Include changes 
from CHANGELOG.rst.
- [ ] Create a pypi package from the master branch. 
      (Use `pip install build` and `python -m build .`)
- [ ] Install the build packages to see if they work.
- [ ] Push tag to remote.
- [ ] Push tested packages to pypi.
- [ ] merge `main` branch back into `develop`.
- [ ] Create a new release on GitHub.
- [ ] Update the package on bioconda.
