#!/bin/sh  
# \
exec tclsh "$0" "$@"

set dir [pwd]

set version 4.10.0
set name spiffe-$version
puts "Building $name RPM"

exec ./rpmdev-setuptree
exec cp -f spiffe.spec $env(HOME)/rpmbuild/SPECS/
exec rm -rf $env(HOME)/rpmbuild/BUILD/$name
exec mkdir $env(HOME)/rpmbuild/BUILD/$name
set binFiles [glob ../bin/Linux-x86_64/*]
foreach f $binFiles {
  exec chmod a+rx $f
  exec chmod a-w $f
  exec cp -f $f $env(HOME)/rpmbuild/BUILD/${name}/
}
cd $env(HOME)/rpmbuild/BUILD
exec tar -cvf ../SOURCES/${name}.tar $name
exec rm -f ../SOURCES/${name}.tar.gz
exec gzip -9 ../SOURCES/${name}.tar
cd ../SPECS
if {[catch {exec rpmbuild -bb --quiet --clean --target x86_64 \
                 --buildroot $env(HOME)/rpmbuild/BUILDROOT spiffe.spec} results]} {
}
exec rm -f ../SOURCES/${name}.tar.gz
puts $results


