# 446-simulate-anneal
Global optimization solution for ENEE446 Assignment 1

Run 'make' on a Mac or Linux terminal to build and './adj' to run

Some smart people came up with a global optimization algorithm that accounts for the pitfalls of the greedy method called simulated annealing. In SA, the algorithm starts with a temperature that cools, and a probability acceptance function that takes in the energy (# of unique pairs) of the neighbors and current states. Unlike the greedy algorithm, which always picks the best neighbor, the SA algorithm weights the best neighbor by the probability acceptance function, so that even worse neighbors may be chosen in an effort to find an even better peak.

http://en.wikipedia.org/wiki/Simulated_annealing

I was able to achieve pretty good results with this alg:

Simulated Annealing method: 197 unique pairs @ 1M samples
Greedy method: 182 unique pairs @ 1M samples
Random sampling: 154 unique pairs @ 1M samples
