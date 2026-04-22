# LOOP-BUFFER-7.0
PROBLEM STATEMENT-
The project aims to develop a computational tool that identifies genetic disease risks within a population by comparing individual DNA sequences against a statistically derived Reference Genome and specific Disease Markers. The system must calculate health severity scores by weighing sequence similarity against the presence of known DNA mutations. 

INTRODUCTION-
The application is a Java-based genomic scanner designed to automate the detection of disorders. It processes large-scale DNA datasets to determine how closely an individual's genetic makeup aligns with a healthy baseline and flags specific sequences that match known disease causing markers. It provides a visual dashboard for risk assessment and clinical monitoring, including data of people who are infected at different severity levels.
•	Reference Genome: A consensus sequence representing the most frequent genetic pattern found in a healthy population.
•	Disease Markers: Specific sub-sequences located at fixed positions in the genome that indicate the presence of a disorder.
•	Genetic Mutation: In this context, a single nucleotide- dna base, mismatch within a disease marker sequence.
•	Risk Severity: A metric calculated by combining the overall genetic similarity to a healthy reference and the total count of detected diseases
TECH STACK- 
Language: Java
GUI Framework: Java Swing 
Data Format: CSV (for healthy datasets, population member datasets, and disease markers)

DATA STRUCTURES AND ALGORITHMS-
1.	Sliding Window:   To know what a "normal" DNA sequence looks like, the program first builds a Reference Genome using k-mer technique The program looks at small "words" (k-mers) of 4 DNA bases (like "ATGC") at every single position across all healthy samples. At position #1, it counts which 4-base "word" appears most often across 300 people. If "ATGC" is the most common at that spot, it becomes part of the reference. Instead of comparing long strings of text, it converts these 4-base words into numbers (integers). Comparing numbers is much faster for a computer than comparing text characters.

2.	Dynamic Programming Table-   Once the "Healthy Baseline" is ready, the program compares a new person's DNA to it to see how similar they are. This uses a famous algorithm called Longest Common Subsequence (LCS). DNA is like a long sentence. Over time, some letters might be deleted or added. LCS finds the longest string of letters that appear in both the healthy reference and the person's DNA in the same relative order, even if they aren't perfectly side-by-side. The program uses a DP Table (a 2D grid) to solve this. It breaks the big DNA string into tiny pieces, solves the similarity for those pieces, and stores the results in the table so it doesn't have to re-calculate them. If the LCS is 450 letters long and the total sequence is 500, the person is 90% similar to the healthy reference.


3.	Hamming Distance Logic-   The program also looks for specific Disease Markers. Each marker has a specific sequence and a specific location where it is known to appear. The program jumps to the exact position (e.g., position 120) and checks if the DNA there perfectly matches the disease sequence (e.g., "TTCG"). Even if it's not a perfect match, the program counts the mismatches.
0 Mismatches: The disease is definitely present (Exact Match).
1 Mismatch: This is flagged as a "Mutation" - the disease might be developing or present in a variant form.

4.	Weighted Algorithm-   It uses a Weighted Algorithm to decide how worried a doctor should be. It assigns a pointer (score) to the similarity % and another pointer to the number of diseases found. It gives more importance (70% weight) to the specific diseases found and less importance (30% weight) to the general similarity. This ensures that even if someone is 99% "normal," finding one dangerous disease marker will still trigger a high-risk warning.
