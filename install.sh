#!/bin/sh
make
mkdir -p ./bin/bin/
cp ./thirdparty/* ./bin/bin/
chmod +x ./bin/bin/*
chmod +x ./bin/kmergo
