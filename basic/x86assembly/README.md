# Basic comparsion of C code with it's corresponding assembly code

## Compliation 

Compiled using `gcc`

```console
$ gcc --version
Apple clang version 14.0.0 (clang-1400.0.29.202)
Target: x86_64-apple-darwin21.6.0
Thread model: posix
```

## Quick C -> Assembly -> Executable -> Final Execution

```bash
./C2AssemblyAndRunExecutable.sh <NAME OF C FILE>
```

```console
$ ./C2AssemblyAndRunExecutable.sh helloworld.c
hello, world!
```

## Quick C -> Assembly

```bash
./C2Assembly.sh <NAME OF C FILE>
```

```console
$ ./C2Assembly.sh helloworld.c
$ cat ./helloworld.s
```

## Quick Assembly -> Executable -> Final Execution

```bash
./AssemblyRunExecutable.sh <NAME OF Assembly FILE>
```

```console
$ ./AssemblyRunExecutable.sh helloworld.s
hello, world!
```