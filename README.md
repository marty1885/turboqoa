# TurboQOA

From scrtch, realtime [Quite Ok Audio](https://qoaformat.org/) codec with streaming support written in pure C. I wrote this thing because I don't want to interdace with the complexity of OGG/OPUS and I don't care that much about audio quanlity and/or compression ratio. QOA's 5x compression is good enough.

## Features

* QOA de-/encoder
* **Streaming** support for both de-/encoder
* Soft realtime algorithm design (could be hard realtime if setup right)
* Dependency free

## TODO:

- [x] Flush content on close of encoding
