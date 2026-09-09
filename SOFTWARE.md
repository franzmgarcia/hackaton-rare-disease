# Software

Python 3.14.0. Original execution: Darwin arm64.

huggingface_hub==1.30.0  
numpy==2.5.3  
scipy==1.18.1  
matplotlib==3.11.1  
pysam==0.24.1  
mappy==2.31  
whatshap==2.8  
pyahocorasick==2.3.1  

External: bcftools 1.22, C++17 compiler and zlib. Build mappy from source if a wheel has unresolved library linkage. No measured end-to-end benchmark. Package generation additionally used reportlab 5.0.1 and pypdf 6.17.0.
