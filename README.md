# QOP

QOP<sup>1</sup> is a quantum variational optimizer simulator written in C.

## Getting Started

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

Running `main` isn't too interesting, so you might want to start by reading the header files, where all declarations are documented<sup>2</sup>, and the source files function comments, if you're interested in the implementation details.

---

<sup>1</sup> pronounced qwop<sup>(kwɒp)</sup>, [like the game](http://www.foddy.net/Athletics.html).

<sup>2</sup> the project attempts to follow [Google's style guide](https://google.github.io/styleguide/cppguide.html).
