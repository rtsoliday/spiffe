Summary:	Accelerator code
Name:		spiffe
License:	EPICS Open license http://www.aps.anl.gov/epics/license/open.php
Group:		Applications/Databases
URL:		http://www.aps.anl.gov/asd/oag/oaghome.shtml
Packager:	Robert Soliday <soliday@aps.anl.gov>
Prefix:		%{_bindir}
Autoreq:	0
Version:	4.10.0
Release:	1
Source:		spiffe-4.10.0.tar.gz


%define debug_package %{nil}
%undefine __check_files
%description
Binary package for Spiffe. Spiffe is a fully electromagnetic 2-1/2
dimensional particle-in-cell code for simulation of rf guns
and similar systems with cylindrical symmetry.

%prep
%setup

%build
%install
mkdir -p %{buildroot}%{_bindir}
install -s -m 755 spiffe %{buildroot}%{_bindir}/spiffe

%files

%{_bindir}/spiffe

