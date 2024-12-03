#!/bin/bash
module load fastqc/0.11.9

fastqc */*.fq.gz -o FASTQC
