This package is meant as a simple way to use the imsrg++ code from Ragnar Stroberg (https://github.com/ragnarstroberg/imsrg) on HPC clusters.

It generates scripts to submit on SLRUM systems using the docker image of the python bindings from the imsrg++ code. The docker image is available at https://hub.docker.com/repository/docker/abelley/pyimsrg/. The image first needs to be converted to a singularity image using:

```console
foo@bar:~$ singularity pull pyimsrg.sif docker://abelley/pyimsrg
```

The jobs submission scripts uses this image to call the python bindings for the IMSRG without needing to install the code base which can be tricky to do on a new cluster.

This code reuses part of scripts from Takayuki Miyagi which can be found at https://github.com/Takayuki-Miyagi/mylib_python.