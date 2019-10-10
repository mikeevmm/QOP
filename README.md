# QOP

QOP¹ is a quantum variational optimizer simulator written in C, with a Python interface.

## Getting Started

## C

In order to build the current codebase, just clone the repository and run `make all`:

```bash
git clone git@bitbucket.org:MikeEVMM/qop.git
cd qop
make all
```

This will compile all the source (`src/`) and headers (`include/`) into objects under `obj/`, and link them into an executable (`main.out`).

You can then run the executable (with entry point at `main.c:main()`) by either explicitly calling

```bash
./main.out
```

or

```bash
make run
```

Running `main` isn't too interesting, so you might want to start by reading the header files, where all declarations are documented², and the source files function comments, if you're interested in the implementation details.

## Python

You can also build and install the Python interface library (`qop`) by running

```bash
# Clone the library if you haven't done so already
# git clone https://MikeEVMM@bitbucket.org/MikeEVMM/qop.git
# cd qop
make py_build
```

This will install the library with the `--user` flag set; if you want to disable this flag, set the `no_user` variable:

```bash
no_user=y make py_build
```

Or if you don't want to install the library at all (in which case it can be found under `build/`), set the `no_install` variable:

```bash
no_install=y make py_build
```

---

¹ pronounced qwop (kwɒp), [like the game](http://www.foddy.net/Athletics.html).

² the project attempts to follow [Google's style guide](https://google.github.io/styleguide/cppguide.html).
