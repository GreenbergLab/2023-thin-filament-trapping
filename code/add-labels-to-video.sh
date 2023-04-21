#!/usr/bin/env nix-shell
#!nix-shell -I nixpkgs=https://github.com/NixOS/nixpkgs/archive/nixos-22.11.tar.gz
#! nix-shell -i bash -p ffmpeg_5-full

ffmpeg -i ../data/motility/apng-videos/biotin-actin-pCa-9.mkv -vf "drawtext=expansion=none:font='Helvetica':text='RTF with 10% biotinylated actin - pCa 9':fontcolor=white:fontsize=30:box=1:boxcolor='#ff4d4dff':boxborderw=10:x=10:y=10" -codec:a copy ../data/motility/apng-videos/biotin-actin-pCa-9-labeled.mkv
