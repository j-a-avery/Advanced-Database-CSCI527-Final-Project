# K-Means
Final project for CSCI-527 Advanced Database.

This readme is a reworked version of the final report submitted with the project.

## 1. Contribution of each member

Due to small class size and precautions taken at the end of the semester due to the COVID-19 pandemic, I was a one-person group.

## 2. Brief implementation summary

The project was implemented in C++ 17. It would probably work compiled as C++ 14 or C++ 11, however several features of modern C++ were used (in particular, functional programming and data structures from the standard library), so it will not compile under C++ 98.

Initially it was planned to implement a matrix template class with a very high-level matrix API, however, this approach expanded the scope of the project far beyond what was required, and introduced a number of nasty memory bugs that would not have been practical to solve in the time alotted. Instead, a more naive, "school-exercise" approach, with more limited abstraction, was deemed to be more practical.

### Distance metrics 

Distance metrics have become a sort of accidental specialty in my degree, so six distance metrics were attempted:

#### Euclidean distance
  
  <img src="https://latex.codecogs.com/png.latex?\textrm{Euclidean-distance}(D,&space;C)&space;=&space;\left\|D&space;-&space;C\right\|_2&space;=&space;\sqrt{\sum_{d&space;\in&space;D,&space;c&space;\in&space;C}&space;\left(d&space;-&space;c\right)^2&space;}" title="\textrm{Euclidean-distance}(D, C) = \left\|D - C\right\|_2 = \sqrt{\sum_{d \in D, c \in C} \left(d - c\right)^2 }" />
  
#### Manhattan distance
  
  <img src="https://latex.codecogs.com/png.latex?\textrm{Manhattan-distance}(D,&space;C)&space;=&space;\sum_{d&space;\in&space;D,&space;c&space;\in&space;C}\left|&space;d&space;-&space;c&space;\right|" title="\textrm{Manhattan-distance}(D, C) = \sum_{d \in D, c \in C}\left| d - c \right|" />

#### Chebyshev distance

  <img src="https://latex.codecogs.com/png.latex?\textrm{Chebyshev-distance}(D,&space;C)&space;=&space;\max_{d&space;\in&space;D,&space;c&space;\in&space;C}\left|&space;d&space;-&space;c&space;\right|" title="\textrm{Chebyshev-distance}(D, C) = \max_{d \in D, c \in C}\left| d - c \right|" />

#### Angular cosine distance

  Cosine similarity is defined as 

  <img src="https://latex.codecogs.com/png.latex?\textrm{cosine-similarity}(D,&space;C)&space;=&space;\frac{D&space;\cdot&space;C}{\left\|&space;D&space;\right\|_2&space;\left\|&space;C&space;\right\|_2}" title="\textrm{cosine-similarity}(D, C) = \frac{D \cdot C}{\left\| D \right\|_2 \left\| C \right\|_2}" />

  that is, the dot product of D and C, divided by the product of the L2 norms of D and C. This can be converted to distance with 

  <img src="https://latex.codecogs.com/png.latex?\textrm{cosine-distance}(D,&space;C)&space;=&space;\frac{^{1}/_{\textrm{cosine-similarity}\left(D,&space;C\right&space;)}}{\pi}&space;\\" title="\textrm{cosine-distance}(D, C) = \frac{^{1}/_{\textrm{cosine-similarity}\left(D, C\right )}}{\pi} \\" />

#### Correlation coefficient as distance

  [1] suggests that correlation coefficient can be used as a distance measure for applications such as *k*-means. The simplest statement of the correlation coefficient formula is

  <img src="https://latex.codecogs.com/png.latex?\rho&space;=&space;\textrm{cosine-similarity}\left(\widetilde{D},&space;\widetilde{C}\right)" title="\rho = \textrm{cosine-similarity}\left(\widetilde{D}, \widetilde{C}\right)" />

  where <img src="https://latex.codecogs.com/png.latex?\inline&space;\tilde{X}&space;=&space;X&space;-&space;\bar{X}\mathbf{1}" title="\tilde{X} = X - \bar{X}\mathbf{1}" />. Since correlation coefficient values fall in the range (-1, 1), they can easily be converted to a distance metric with <img src="https://latex.codecogs.com/png.latex?\inline&space;1&space;-&space;\left(\frac{\rho}{2}&space;-&space;\frac{1}{2}\right)" title="1 - \left(\frac{\rho}{2} - \frac{1}{2}\right)" />.

#### Real-valued Jaccard distance

  Jaccard similarity is typically presented as a measurement of similarity between two sets:

  <img src="https://latex.codecogs.com/png.latex?\frac{\left|S&space;\cap&space;R\right|}{\left|&space;S&space;\cup&space;R&space;\right|}&space;=&space;\frac{\left|&space;S&space;\cap&space;R&space;\right|}{\left|S\right|&space;&plus;&space;\left|R\right|&space;-&space;\left|S&space;\cap&space;R\right|}" title="\frac{\left|S \cap R\right|}{\left| S \cup R \right|} = \frac{\left| S \cap R \right|}{\left|S\right| + \left|R\right| - \left|S \cap R\right|}" />

  However, [2] suggests a form generalized to vectors, which can be used as a similarity metric:

  <img src="https://latex.codecogs.com/png.latex?\textrm{real-valued-Jaccard-similarity}\left(D,&space;C&space;\right)&space;=&space;\frac{D&space;\cdot&space;C}{D^2&space;&plus;&space;C^2&space;-&space;D\cdot&space;C}" title="\textrm{real-valued-Jaccard-similarity}\left(D, C \right) = \frac{D \cdot C}{D^2 + C^2 - D\cdot C}" />

  where of course <img src="https://latex.codecogs.com/png.latex?\inline&space;X^2&space;=&space;X&space;\cdot&space;X" title="X^2 = X \cdot X" />.

Two of these measures &mdash; angular cosine distance and correlation coefficient distance &mdash; had memory errors; since four distance metrics were more than enough to meet the requirements of the project, these two were excluded from the final analysis.

### Centroid selection schemes

In addition to the six distance metrics, two centroid selection schemes were used:
    
1. Select centroids at random from among the existing data points
2. Generate centroids at random from within the range of the data

### Demonstration

To demonstrate the multiple distance measurements, a number of modern C++ features, including type aliasing, functional programming capabilities, and data structures from the standard library were used. Since each of the measurement functions has the same type signature &mdash; `matrix& f(matrix&, matrix&)` &mdash; that type signature was defined as a type `measure_function`, with `pmeasure_function` used for pointers to `measure_function`s. Similarly, since the functions to generate centroids have type `matrix& f(int, matrix&)`, type aliases were defined for this signature and its pointer. Finally, tuples of a function and its name were defined using `std::pair`: `std::pair<std:;string, pmeasure_func> named_measure_function` and `std::pair<std:;string, pcentroid_generator> named_centroid_generator`. Since the name and function are accessed with `.first` and `.second` respectively, wrapper functions were created to extract the name and function from the pair. (An alternative approach which would have been more idiomatic for older versions of C++ would have been to wrap the name and function pointer in a `NamedFunction` class, rather than a `std::pair`; the functional approach was much more straightforward and much simpler to implement with the deadline approaching.)

For testing and demonstration, two `std::map`s were created, one populated with the distance measures, and the other populated with the centroid generators. These were then looped through for values of *k* between 3 and 5, inclusive, each of the distance measures, and each centroid policy, with 5 runs for each.


## Post-Mortem when Posting to GitHub

From a programming language standpoint, I feel somewhat like an immigrants' child who has forgotten his parents' tongue (C++) in favor of the native language of his new homeland (Python). I've used C++ in CS1 and data structures courses at every step of my education, but never further than that. I've done projects much larger than this one in Python, Java, and C#, but prior to this project my C++ experience was limited to smaller programming exercises.

K-means is a fairly simple algorithm that can be found in just about any machine learning or data mining textbook, and even some linear algebra books, so implementing it in Python can be done in an afternoonor less, but this project specifically required C++ without external libraries. (Even using the C++17 standard library felt transgressive, as if I was following the letter of the law but not its spirit, but using FP principles and not reinventing certain wheels streamlined the project considerably.) But, of course, what takes hours in Python (or, quite frankly, any other language) takes weeks or longer in C++, especially when the project requires not just the wheel to be reinvented, but the ax and the saw as well. (This is a pattern I've noticed before &mdash; a project that took two months in C++, despite the project *explicitly calling for code to be copied from the textbook*, took four days to redo in Go, from scratch, including the time required to learn the language.)

What I originally wanted to do for this project was to implement a very high-level, almost Pythonic, *e.g.* the code to calculate Euclidean distances would look something like the following:

```c++
template <typename Number>
matrix<double>& Euclidean_distances(matrix<Number>& D, matrix<Number>& C) {
	return (D.as_double() - C.as_double).norm(DIM::ROWS);
}
```

Parts of this approach worked very well; other parts caused memory errors. Some of those memory errors were simple to fix, and others caused kernel panics. A week before the deadline I started over. Lesson learned: start simple and work up to more complexity.


## References

[1] S. Boyd and L. Vandenberghe, *Introduction to Applied Linear Algebra: Vectors, Matrices, and Least Squares*. Cambridge University Press, 2018.

[2] C. C. Aggarwal, *Machine Learning for Text*. Springer International Publishing, 2018.
