# Documentation build instructions

To generate the LaTeX documentation in this directory, ensure the required
TeX and helper utilities are installed:

```
apt-get update
apt-get install -y --no-install-recommends \
  texlive-base texlive-latex-base texlive-latex-extra \
  texlive-fonts-recommended texlive-plain-generic \
  texlive-extra-utils tex4ht ghostscript latex2html \
  texlive-font-utils dvipng
```

If an error related to `ca-certificates-java` appears during installation,
remove the OpenJDK packages before retrying:

```
apt-get remove --purge -y openjdk-*-jre-headless default-jre-headless default-jre
```

After installing dependencies, run `make` in this directory to build the docs.

