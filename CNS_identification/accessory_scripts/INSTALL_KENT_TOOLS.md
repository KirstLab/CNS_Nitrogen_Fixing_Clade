# Kent MAF utilities (mafExtract, mafIndex)

`CNS_identification/Snakefile` calls `./mafExtract` and `./mafIndex` from its working directory. These binaries are not committed to the repository because the upstream copies are Linux x86_64 only and add roughly 10 MB to a clone.

## Linux x86_64

Fetch the prebuilt binaries from the CESAR2.0 distribution (the same source already referenced in `CNS_identification/config.yaml`):

```bash
cd CNS_identification
curl -fsSL https://github.com/hillerlab/CESAR2.0/raw/master/tools/mafExtract -o mafExtract
curl -fsSL https://github.com/hillerlab/CESAR2.0/raw/master/tools/mafIndex   -o mafIndex
chmod +x mafExtract mafIndex
```

## macOS or ARM hosts

The CESAR2.0 binaries will not run. Build the equivalent Kent utilities from source from the UCSC Genome Browser source tree (`https://github.com/ucscGenomeBrowser/kent`) and place the resulting `mafExtract` and `mafIndex` executables in `CNS_identification/`.
