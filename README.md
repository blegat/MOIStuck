# Reproduce

Use Julia v1.0, v1.1 or Julia master.
Fun `bug.sh`. If you get:
```
Calling
bug.sh: line 1: 11366 Killed                  timeout -s 9 10 /home/blegat/git/julia-master/julia --color=yes bug.jl
```
It means you hit the bug. If you get instead:
```
Calling
Called
...
```
It means that you simplified it too much.
