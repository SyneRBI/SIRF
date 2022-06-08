## Release Checklist
Set version variable of the release for further steps, e.g. using the format
```
VER=3.2.0-rc.1
```

1. Submodules (within SIRF)
  - [ ] tag data
    + `git submodule update` # **should not complain**
    + `cd data && git tag -a v$VER -m "version $VER" && git push --tags && cd ..`
    + `git status`  # **should not list `data` as modified**
2. SIRF
  - [ ] update `CHANGES.md`
  - [ ] update `NOTICE.txt`
  - [ ] update `CITATION.cff`
  - [ ] update version numbers in [SIRF/CMakeLists.txt](https://github.com/SyneRBI/SIRF/blob/master/CMakeLists.txt)
  - [ ] update version numbers in the [doc/UsersGuide.md](https://github.com/SyneRBI/SIRF/blob/master/doc/UserGuide.md) etc
  - [ ] run all demos
  - [ ] run ctest
  - [ ] `git push`
  - [ ] check Travis
  - [ ] run doxygen, check, send files for uploading and update the doxygen link on Wiki
  - [ ] `git tag -a v$VER -m "version $VER"`
  - [ ] `git push origin v$VER`
  - [ ] if final release only: create release on https://github.com/SyneRBI/SIRF/releases/new linked to the tag and copying all CHANGES.md since last release (as listed on https://github.com/SyneRBI/SIRF/releases)
3. SuperBuild
  - [ ] update `CHANGES.md`
  - [ ] update `NOTICE.txt`
  - [ ] update `CITATION.cff`
  - [ ] update `SIRF-Superbuild/version_config.cmake` with new `SIRF_TAG` (and `STIR_TAG` etc if necessary)
  - [ ] update version number in [VM_version.txt](https://github.com/SyneRBI/SyneRBI_VM/blob/master/VM_version.txt)
  - [ ] update `vb.name` in [vagrant/vagrantfile](https://github.com/SyneRBI/SyneRBI_VM/blob/master/vagrant/Vagrantfile)
  - [ ] `git push`
  - [ ] check Travis
  - [ ] `git tag -a v$VER -m "version $VER"`
  - [ ] `git push origin v$VER`
  - [ ] if final release only: create release on https://github.com/SyneRBI/SIRF-SuperBuild/releases/new linked to the tag and copying all CHANGES.md since last release (as listed on https://github.com/SyneRBI/SIRF-SuperBuild/releases)
4. SIRF-Contribs
  - [ ] `git tag -a v$VER -m "version $VER"`
  - [ ] `git push origin v$VER`
5. Virtual Machine
  - [ ] `vagrant up`
  - [ ] Virtualbox Guest Additions
  - [ ] run `first_run.sh` script (gnome settings and zero-fill trick)
  - [ ] check that the serial port is [deselected](https://github.com/SyneRBI/SyneRBI_VM/blob/master/vagrant/README.md#notes-about-ubuntu-box-for-version-100).
  - [ ] export the VM
  - [ ] ctest on VM
  - [ ] run all exercises (download data first)
  - [ ] upload to Zenodo
        * use name like `SIRF 2.2.0.ova` (with space) as Zenodo uses alphabetical ordering probably, add `-release-candidate-1` if necessary (as most people will not know what `rc1` means).
6. SIRF-Exercises (already checked in VM)
  - [ ] `git tag -a v$VER -m "version $VER"`
  - [ ] `git push origin v$VER`
7. Website (if final release)
  - [ ] update Software page (version info, VM etc)
  - [ ] upload doxygen
  - [ ] update link for doxygen in [Wiki](https://github.com/SyneRBI/SIRF/wiki/Software-Documentation)
  - [ ] add news flash
8. Announce
  - [ ] Send email to SyneRBI-DEVEL@JISCMAIL.AC.UK; SyneRBI-USERS@JISCMAIL.AC.UK; add SyneRBI@JISCMAIL.AC.UK for final release
9. Tag wiki
  - [ ] `git clone https://github.com/SyneRBI/SIRF.wiki.git; cd SIRF.wiki` (or pull)
  - [ ] `git tag -a v$VER -m "version $VER"`
  - [ ] `git push origin v$VER`
 