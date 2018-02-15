## Release Checklist
Set version variable of using the format
````
VER=1.0.0-rc.1
```

1. SIRF
 - [ ] update `CHANGES.md`
 - [ ] update `NOTICE.txt`
 - [ ] update version numbers in [SIRF/CMakeLists.txt](CMakeLists.txt)
 - [ ] update version numbers in the [doc/UsersGuide.md](doc/UserGuide.md) etc
 - [ ] `git push` 
 - [ ] check Travis
 - [ ] run doxygen, check and send files for uploading
 - [ ] `git tag -a v$VER -m "version $VER"`
 - [ ] `git push origin v$VER`
 
2. SuperBuild
 - [ ] update `CHANGES.md`
 - [ ] update `NOTICE.txt`
 - [ ] update `SIRF-Superbuild/version_config.cmake` with new `SIRF_TAG`
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
 - [ ] remove screen lock
 - [ ] add English UK keyboard
 - [ ] run the zerofill trick
 - [ ] Virtualbox Guest Additions
 - [ ] export the VM and upload to website
 - [ ] `git tag -a v$VER -m "version $VER"`
 - [ ] `git push origin v$VER`

4. Website
 - [ ] update Software page (version info, VM etc)
 - [ ] upload doxygen
 - [ ] update link in [Wiki](https://github.com/CCPPETMR/SIRF/wiki/Documentation)
 - [ ] add news flash

5. Announce
 - [ ] Send email


 #### zerofill trick
 
 ```
sudo dd if=/dev/zero of=/EMPTY bs=1M
sudo rm -f /EMPTY
```


