## How to use the code

This script is built to be fully interactive so you can test it on any dataset without having to hardcode variables every single time. Just hit the **Run All** button on the cells below and follow the terminal prompts:

1. **Choose your input mode:**
   * Type `y` if you want to manually drop in your own test sequences.
   * Type `n` to just run the default sequence block from a youtube video (Check Reference).
2. **Enter sequences (if you chose `y`):**
   * Type each DNA sequence one by one and hit `ENTER`.
   * When you are finished adding sequences, literally just type `DONE` to lock them all in.
3. **Set the Motif Length ($l$):**
   * The code will ask if you have a specific "Target Motif" you're hunting for (e.g., `GATTACA`). If you type it in, the code automatically sets $l$ to the length of that word and does a quick flex by scoring how closely its final answer matches your target.
   * If you don't have a target and just want to blindly search for a pattern of a certain size, hit `ENTER` to skip, and then type the integer length (e.g., `6`).

**What to expect in the output:**
Because `show_steps=True` is ON, the code is gonna dump a huge **Step-by-Step Breakdown** during its very first iteration. It includes the raw count matrix, the pseudocount probability derivation, and the literal mathematical substitution for the $R(w)$ weights so that we can crosscheck. After this, it spits out the final Log-Odds Matrix and the best aligned motifs across all sequences. (For CrossChecking purposes)

### References & Acknowledgements

* **Primary Algorithm Source:** Rouchka, E. C. (2008). *A Brief Overview of Gibbs Sampling*. [ResearchGate Link](https://www.researchgate.net/publication/2336786_A_Brief_Overview_of_Gibbs_Sampling)
* **Biological Context & Framework:** Oregon State University. *Applied Bioinformatics - Chapter 2: Sequence Motifs*. [OER Link](https://open.oregonstate.education/appliedbioinformatics/chapter/chapter-2-sequence-motifs/)
* **Youtube:** [Youtube Link](https://youtu.be/KqquwB_kHvU)

### Algorithm Considerations & Drawbacks
Based on the lecture notes and the actual way this code runs, there are some tradeoffs we have to accept:
* **The Good:** The core algorithm runs in O(N) linear run-time. Because of the stochastic roulette-wheel sampling, it actively avoids the local optimum problem. We can also view suboptimal results if we want to trace the iterations.
* **The Bad (Algorithm Limits):** It requires relatively large sets of sequences to work properly. Even with our probabilistic sampling, it is not guaranteed to reach an optimal value. Also, we must specify an estimate for the motif width upfront. If the actual biological motif is 12 bases long and I tell the code to look for an $l$-mer of 6, it completely misses the big picture.
* **The Drawbacks of My Specific Code:** Because I wrote this using native Python dictionaries and nested loops to make the step-by-step math printable for the manual trace, it's pretty unoptimized for speed. If we threw a massive real-world genomic dataset at this, it would lag terribly compared to a vectorized NumPy/C++ version. Also, to fight the fact that Gibbs isn't actually guaranteed to reach the optimal value, I brute-forced it by adding a massive `num_restarts` outer loop, which multiplies the run-time heavily just to secure a safe grade.

### Deep Learning Crossovers (Why the code is structured this way)
While this is a bioinformatics algorithm, I actually ported over several concepts directly from deep learning to optimize how the Gibbs Sampler converges:

1. **Early Stopping / Patience (`patience_counter`):** Instead of looping blindly for $N$ iterations, I implemented an early stopping callback. If the model's convergence score ($L$) fails to improve after 50 continuous epochs (iterations), training halts. This prevents overfitting to a local motif and saves compute time.
2. **Laplace Smoothing (Pseudocounts):** I added `+1.0` to the background and profile matrices. In NLP/Deep Learning embeddings, an unseen token will crash the model with a zero-probability multiplication. Slapping a pseudocount here ensures our mathematical chain doesn't break if a sequence is missing a specific nucleotide.
3. **Log-Likelihood (`math.log2`):** Multiplying raw probabilities over and over causes severe floating-point underflow (eventually just rounding to 0.0 in Python). Using Log-Odds allows us to *add* scores instead, identical to how neural networks use Cross-Entropy or LogSoftmax loss functions for stability.
4. **Softmax-style Sampling (Roulette Wheel):** Gibbs sampling gets trapped easily if you just use `argmax` (greedy decoding) to pick the best l-mer. I normalized the weights and sampled probabilistically—essentially the same logic as using a temperature parameter in Softmax to allow the model to "explore" alternative alignments.
5. **Random Weight Initialization (`num_restarts`):** A bad random start means getting stuck in a local minimum. Just like re-initializing weights in Stochastic Gradient Descent, running multiple random restarts guarantees the algorithm eventually finds the global optimal motif.

### Integration with Previous Bioinformatics Coursework
To tie this into the past assignments, the core of this algorithm builds directly on foundational bioinformatics concepts taught in class so far:
1. **Position Weight Matrices (PWMs / PSSMs):** The `create_profile` function is essentially just building a PWM. We tally the raw frequencies of the motif, apply pseudocounts, and normalize them into probabilities to score unseen sequences, exactly like we did in the basic motif evaluation labs.
2. **Gapless Local Multiple Sequence Alignment (MSA):** Unlike exact dynamic programming methods (like Needleman-Wunsch or Smith-Waterman) used for pairwise alignment, this Gibbs sampler performs a heuristic, stochastic *local* alignment across multiple sequences simultaneously, strictly without gaps.
3. **Zero-Order Background Model:** The `get_background` function calculates a global probability for A, C, G, and T. This assumes each nucleotide's appearance is independent of its neighbors, serving as our zero-order Markov null model to test the motif's statistical significance against.



> **Note:** AI (Google Gemini) was used as a coding assistant to translate the course's mathematical formulas and pseudocode into Python. All core algorithm logic is derived strictly from the referenced course materials.
