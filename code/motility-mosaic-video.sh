#!/usr/bin/env nix-shell
#!nix-shell -I nixpkgs=https://github.com/NixOS/nixpkgs/archive/nixos-22.11.tar.gz
#! nix-shell -i bash -p ffmpeg_5-full
ffmpeg \
	-i ../data/motility/apng-videos/mosaic/regular-actin-pCa-9-labeled-cropped.mkv \
	-i ../data/motility/apng-videos/mosaic/regular-actin-pCa-625-labeled-cropped.mkv \
	-i ../data/motility/apng-videos/mosaic/regular-actin-pCa-4-labeled-cropped.mkv \
	-i ../data/motility/apng-videos/mosaic/biotin-actin-pCa-9-labeled-cropped.mkv \
	-i ../data/motility/apng-videos/mosaic/biotin-actin-pCa-625-labeled-cropped.mkv \
	-i ../data/motility/apng-videos/mosaic/biotin-actin-pCa-4-labeled-cropped-scale-bar.mkv \
	-filter_complex "nullsrc=size=2724x1352 [base];\
    [0:v] setpts=PTS-STARTPTS, scale=908x676 [upperleft];\
    [1:v] setpts=PTS-STARTPTS, scale=908x676 [uppermiddle];\
    [2:v] setpts=PTS-STARTPTS, scale=908x676 [upperright];\
    [3:v] setpts=PTS-STARTPTS, scale=908x676 [lowerleft];\
    [4:v] setpts=PTS-STARTPTS, scale=908x676 [lowermiddle];\
    [5:v] setpts=PTS-STARTPTS, scale=908x676 [lowerright];\
    [base][upperleft] overlay=shortest=1 [tmp1];\
    [tmp1][uppermiddle] overlay=shortest=1:x=908 [tmp2];\
    [tmp2][upperright] overlay=shortest=1:x=1816 [tmp3];\
    [tmp3][lowerleft] overlay=shortest=1:x=0:y=676 [tmp4];\
    [tmp4][lowermiddle] overlay=shortest=1:x=908:y=676 [tmp5];\
    [tmp5][lowerright] overlay=shortest=1:x=1816:y=676" \
	-c:v libx264 ../data/motility/apng-videos/mosaic/motility-mosaic.mkv
