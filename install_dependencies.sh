#!/bin/bash

if ! command -v cellranger >/dev/null 2>&1; then
    mkdir -p resources

    if [ ! -f "resources/cellranger-9.0.1.tar.gz" ]; then
        wget -O resources/cellranger-9.0.1.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-9.0.1.tar.gz?Expires=1750008677&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=YKDCWaPOAfekE26khc1L8S2YrvikKB1u8SE6tBYq3QFGEE7jroiOTLxPNb6sMhKVG8uhVWUTcUyky4VwH3aVwXKF76bw3wwOJBHH3Tb8z4nW~Y9wCc4tmd0QVLFrePtB2cycgw8LFt~p~ZSc-zlK--TDycchZrJo0AzeeYRuoNapj6BukUPa659-qdytHP1YyJLhLDKBAdeKL3jn9iJUOCHyiuZRP1lOx1Iu3bNd0XltvOWT~ZFzqKKug8epQgcvbwGV3nESwKBjwRAhXvMfOajpUDYNFzztpLFEDcaZgKk34Rh0Yxv1~PmJFiH0Bf78OHSv8cljIENZVEUkAbQBzw__"
    fi

    if [ ! -d "resources/cellranger-9.0.1" ]; then
        tar -xzvf resources/cellranger-9.0.1.tar.gz -C resources
    fi

    if ! grep -q 'resources/cellranger-9.0.1' ~/.zshrc; then
        echo 'export PATH="$PWD/resources/cellranger-9.0.1:$PATH"' >> ~/.zshrc
        export PATH="$PWD/resources/cellranger-9.0.1:$PATH"
    fi

    echo "cellranger 9.0.1 installed and configured in resources/."
else
    echo "cellranger is already installed and in PATH."
fi

if [ -f resources/Homo_sapiens.GRCh38.113.chr.gtf.gz ] && [ ! -f resources/Homo_sapiens.GRCh38.113.chr.gtf ]; then
    gunzip -c resources/Homo_sapiens.GRCh38.113.chr.gtf.gz > resources/Homo_sapiens.GRCh38.113.chr.gtf
fi

if [ ! -f resources/GRCh38.fa ]; then
    wget -O resources/GRCh38.fa.gz "https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
    gunzip resources/GRCh38.fa.gz
fi

cd resources

cellranger mkref --genome=GRCh38_113 \
  --fasta=GRCh38.fa \
  --genes=Homo_sapiens.GRCh38.113.chr.gtf \
  --nthreads=32

cd ..
