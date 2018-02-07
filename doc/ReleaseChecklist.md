## Release Checklist

1. SIRF
 - [ ] update `CHANGES.md`
 - [ ] update version numbers in `SIRF/CMakeLists.txt`
 - [ ] run doxygen and upload files
 
2. SuperBuild
 - [ ] update `CHANGES.md`
 - [ ] update `SIRF-Superbuild/version_config.cmake`
 - [ ] run tests (and check Travis)

3. Virtual Machine

 - [ ] run `vagrant up`
 - [ ] run the zerofill trick
 - [ ] export the VM and upload to website
 
 #### zerofill trick
 
 ```
sudo dd if=/dev/zero of=/EMPTY bs=1M
sudo rm -f /EMPTY
```

### Documentation
