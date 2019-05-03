## Release Checklist
Set version variable of the release for further steps, e.g. using the format
```
VER=1.0.0-rc.1
```

1. SIRF
 - [ ] update `CHANGES.md`
 - [ ] update `NOTICE.txt`
 - [ ] update version numbers in [SIRF/CMakeLists.txt](https://github.com/CCPPETMR/SIRF/blob/master/CMakeLists.txt)
 - [ ] update version numbers in the [doc/UsersGuide.md](https://github.com/CCPPETMR/SIRF/blob/master/doc/UserGuide.md) etc
 - [ ] run all demos
 - [ ] run ctest
 - [ ] `git push` 
 - [ ] check Travis
 - [ ] run doxygen, check, send files for uploading and update the doxygen link on Wiki
 - [ ] `git tag -a v$VER -m "version $VER"`
 - [ ] `git push origin v$VER`
 
2. SuperBuild
 - [ ] update `CHANGES.md`
 - [ ] update `NOTICE.txt`
 - [ ] update `SIRF-Superbuild/version_config.cmake` with new `SIRF_TAG` (and `STIR_TAG` etc if necessary)
 - [ ] `git push`
 - [ ] check Travis
 - [ ] `git tag -a v$VER -m "version $VER"`
 - [ ] `git push origin v$VER`

3. Virtual Machine

 - [ ] update version number in [VM_version.txt](https://github.com/CCPPETMR/CCPPETMR_VM/blob/master/VM_version.txt)
 - [ ] update `vb.name` in [vagrant/vagrantfile](https://github.com/CCPPETMR/CCPPETMR_VM/blob/master/vagrant/Vagrantfile)
 - [ ] update `CHANGES.md`
 - [ ] update `NOTICE.txt`
 - [ ] `git push`
 - [ ] `vagrant up`
 - [ ] Virtualbox Guest Additions
 - [ ] run `first_run.sh` script (gnome settings and zero-fill trick)
 - [ ] [deselect](https://github.com/CCPPETMR/CCPPETMR_VM/blob/master/vagrant/README.md#notes-about-ubuntu-box-for-version-100) the serial port.
 - [ ] export the VM
 - [ ] ctest on VM
 - [ ] upload to website
 - [ ] `git tag -a v$VER -m "version $VER"`
 - [ ] `git push origin v$VER`

4. Website
 - [ ] update Software page (version info, VM etc)
 - [ ] upload doxygen
 - [ ] update link for doxygen in [Wiki](https://github.com/CCPPETMR/SIRF/wiki/Software-Documentation)
 - [ ] add news flash

5. Announce
 - [ ] Send email

6. Tag wikis
 - [ ] `git clone https://github.com/CCPPETMR/SIRF.wiki.git; cd SIRF.wiki` (or pull)
 - [ ] `git tag -a v$VER -m "version $VER"`
 - [ ] `git push origin v$VER`
 - [ ] `git clone https://github.com/CCPPETMR/CCPPETMR_VM.wiki.git; cd CCPPETMR_VM.wiki` (or pull)
 - [ ] `git tag -a v$VER -m "version $VER"`
 - [ ] `git push origin v$VER`

