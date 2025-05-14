
Notes for the book *Deep Learning*: 

Goodfellow, I., Bengio, Y., & Courville, A. (2017). Deep Learning. The MIT Press. 



## Chapter 1: Intro

### ‘Deep’

-   No consensus about how much depth a model requires to qualify as ‘deep’ due to different measuring methods of model depth (longest path of steps in a sequential operation (computational graph) vs. layers of related *concepts* (probabilistic modeling graph))
-   Layers: visible layer - contains the variables that *we* are able to *observe*, eg. the input; hidden layers.
-   AI -> machine learning -> representation learning -> deep learning

![1744852728384.png](https://img.picui.cn/free/2025/04/17/680056fbbe642.png)

-   **Representation learning**: differs from classic ML as representations (features) are learned by the system rather than manually designed. Eg. the **autoencoder** - encoder (input data -> representation) + decoder (the new representation -> back into the original format).

  

-   Deep learning: building complex representations out of simpler ones. Eg. **multilayer perceptron (MLP)**, just a mathematical *function* (formed by composing many *simpler functions*) that maps some set of input values to output values. Each application of a different function provides a new representation of the input.

-   Three waves of NN research:

	-   ‘Cybernetics’: McCulloch-Pitts Neuron (early model of brain function); perceptron (first model that could learn the weights *w*); ADALINE (already applies the algorithm stochastic gradient descent).

		-   **Linear models**, f(x, w) = x<sub>1</sub>w<sub>1</sub> + … + x<sub>n</sub>w<sub>n</sub> -> y. Fails to solve *nonlinear* problems, eg. the XOR function, where f([0,1], w) = 1 and f([1,0], w) = 1 but f([1,1], w) = 0 and f([0,0], w) = 0.

	-   ‘Connectionism’: a large number of simple computational units can achieve intelligent behavior when *networked together*, an idea inspired by the Hebbian Principle (simplified as that ‘fire together wire together’ thingy).

		-   Distributed representation: representing each input by separate features.
		-   Back-propagation.

	-   ‘Deep learning’: leveraging 1. increased dataset sizes (the age of “Big Data”); 2. computational resources to run much larger models.

-   Neuroscience is no longer the predominant guide of DL research. While it did inspire several NN _*_architectures_*,_ it does not offer much guidance for the *learning algorithm* used to *train* these architectures.

	-   Artificial NNs are much more densely connected (more connections *per neuron*) than in many mammalian biological NNs.

  

---  

  

Part 1: Applied Math and ML Basics

> Hoping to grasp a better understanding of these math subjects under the context of DL…

## Chapter 2: Linear Algebra

// Skimmed for most, tho there are always parts i can't understand (some are so very simple, i'm dumb), so i made notes (either still typed in markdown or written) for those to help me,,, ><

  

-   **Arrays…**

	-   Vectors: can be thought of as identifying points in space, with each element giving the *coordinate along a different axis* (eg. for Cartesian systems it’s the x/y- coordinate directed along the i-hat/ j-hat, etc.

		-   The - sign is used to index the complement of a set {}, eg. **x**<sub>-1</sub> is the vector containing all elements of **x**  except for *x<sub>1</sub>*, similar to negative indexing in R. 

	-  Matrices are often viewed as column vectors *stacked* side-by-side, so the columns would specify the 'directions' in a linear combination perspective (elaborated afterwards). The rows can still be regarded as axes in that *space* (then in this case 'space' would mean rowspace; this view got consolidated after reading the last part on PCA...)
		> (i love the word 'space' somehow,,,,, it gives the feeling of extensive possibilities since one can grab virtually any feature they like to be the axis in that *space*...)
		
		- I think it's necessary to point out though, that in data science, datasets (*design matrices*)  often structure *rows* as observations or sample points and columns as features, which is sometimes the transpose of other notions. So ***context matters***!! (just like how the term 'parameter' has differing meanings: in ML it's the model's *internal variables* (i.e. weights *w*, learned from data) not the input variable x; in plain programming it's the formal variable in a function definition)
	- Tensors: an array of numbers with a variable number of axes.
	-   Distinguish the **Hadamard** product of matrices (element-wise!) from  the classic `.matmul()` operation.
	-   **Symmetric** matrices often arise for some function of two arguments that does not depend on the *order* of the arguments: eg. if A  is a matrix of distance measurements, with A<sub>_i,j_</sub>  giving the distance from point _i_ to point _j_, then A<sub>i,j</sub> = A<sub>j,i</sub> because distance functions are symmetric!!
	-   **Orthogonal** matrices have mutually ortho*normal* columns and rows.

  

-   **Ax = b…**

	-   Using the inverse of A to solve the linear equation Ax = b -> x = A^(-1)b is not actually used in practice due to limited precision of representing A^(-1) on a digital computer. Algorithms that make use of the value of b can usually obtain more accurate estimates of x, the vector with unknown variables.
	-   We can think of the *columns* of A  as specifying different directions we can *travel* from the origin  (the point specified by the vector of all zeros), and determine how many ways there are of reaching b. In this view of linear combination, each element of x specifies how far we should travel in *each of these directions*, with x<sub>i</sub>  specifying how far to move in the direction of column i.
	-   For a non-square matrix A (m x n, n > m), it is possible for them to have more than one *set* of m linearly dependent columns (when the rank(A) = m, so that it is possible to find a solution for x with a given b that lies in the column space with m dims). But why does such a matrix have to have *infinitely many* (rather than a *finite* set), if not one, solutions??? (ok i'm dumb about forgetting how this part was taught in class ><,, but reviewing the *nullspace* rlly helps explain this statement and the following part on how the pseudoinverse finds the solution with the minimal L2 norm)
![1744852754168.png](https://img.picui.cn/free/2025/04/17/6800571605a78.png)
	- ***See the Moore-Penrose Pseudoinverse below***.

  

-   **Norms**

	-   Another way to interpret them is *functions that map vectors to non-negative values*.
	-   The **squared** L2 norm is common and convenient to use (simply calculating the dot product of the transpose of x and x, since that's just adding up the squared values of each element… typical dot product behavior), because the derivatives of the *squared* L2 norm with respect to each element of x depend only on the corresponding element of x itself (-> *2x<sub>i</sub>*) , while for the L2 norm (involving the *square root*!!) the derivatives with respect to an element depend on the *entire vector* through the norm in the denominator. (That must also be part of the reason why minimizing the *squared* error is often used..? Since *gradient*-based optimization methods would work more efficiently.)
	-   The squred L2 norm may be undesirable sometimes tho as it increases very slowly near the origin. Thus, in cases when we want to discriminate between elements that are exactly zero and those that are small but still nonzero, the L1 norm (physical distance..) is used, as the L1 norm grows at the same rate (just increases by that distance moved) in all locations.
	-   **Frobenius** norm (denoted by ||A||<sub>F</sub>): just the plain L2 norm but applied to elements of a matrix rather than a vector. An alternative way to write it is sqrt(Tr(AA<sup>T</sup>)), which is very obvious after expanding this expression.

  

-   **Other notations: **

	-   Broadcasting is denoted by the shorthand notation C = A + b.
	-   The ‘turned A’ （upside-down) symbol: denotes ‘For All’.

  

### Matrix Decompositions

**Mathematical objects can be understood better by breaking them into constituent parts, just like decomposing an integer into its prime factors helps us understand that integer’s divisibility.** 

-   Eigendecomposition
	- Eigen- stuff is so ubiquitous in later concepts (PCA, Hessian matrices, etc).  Eigendecomposition (requires square matrices, and for some of them real eigenvalues don't exist) splits complex transformations into simpler scaling operation along *specific directions* that are invariant (these directions are specified by the eigenvectors; the eigenvalues quantify the extent of that scaling. and since the resulting matrix containing the eigenvalues is diagonal (*diagonalization*), this also simplifies computations).
  
	> The eigendecomposition of a real symmetric matrix can also be used to optimize quadratic expressions of the form f(x) = **x<sup>T</sup> Ax** subject to ||x||<sub>2</sub> = 1
	
	- I couldn't understand the meaning of this (after I read the following chapters, it turned out that the terms were only officially introduced then). '**Optimize**' - just means finding the value of the *argument* (i.e., the input variable *x*) that *maximizes or minimizes a function f(x)*, often denoted using an asterisk (eg. f*, x*); `argmax`/`argmin`. '**subject to**' - marks the *constraint* of the optimization (in this case the constraint is the L2 norm of x must be 1, which forces x to lie on the unit sphere). 
	- When A is positive semidefinite, for all x, x<sup>T</sup> Ax >= 0. 
		- I rewatched 3blue1brown's [Eigenvectors and eigenvalues](https://www.youtube.com/watch?v=PFDu9oVAE-g), focusing on the last part about *eigenbasis*, in which it explains how we can consider the transformation by a matrix A in another (usually simpler... since the resulting matrix would be a diagonal matrix) way: changing our grid (i-hat, j-hat) to a new grid where the eigenvectors of A are now the basis vectors. I find it quite helpful to understand our statements here about f(x) = x<sup>T</sup> Ax and why it applies to *all* x for a positive semidefinite matrix A, so again I kept some notes,,,  (The prior 3b1b episode on [change of basis](https://www.youtube.com/watch?v=P2LTAUO1TdA) is also so fun... i'll never forget that back-and-forth translation from 'our language' to 'Jennifer's language', so intuitive lolol)
	
	(I took the note below when reading about the case in which A is a real symmetric matrix (**normal**) to try to understand the corresponding statement from the book; the change-of-basis thingy applies to more general cases of A (non-symmetric))![1744852768585.png](https://img.picui.cn/free/2025/04/17/68005724c1117.png)
	
-   SVD: now applicable to non-square matrices or matrices without real eigenvectors, singular values are pretty much analogous to eigenvalues and singular vectors generalize eigen-concepts. The most useful feature of it is to partially generalize *matrix inversion* to non-square matrices:

	-   **The Moore-Penrose pseudoinverse**

		-   A = UDV<sup>T</sup> (from SVD)  ->  A<sup>+</sup> (the pseudoinverse of A) = VD<sup>+</sup>U<sup>T</sup>
		-  The pseudoinverse of D (D<sup>+</sup>) is obtained by taking the reciprocal of its non-zero elements and then transposing it. (The inverse of a square diagonal matrix is just taking the reciprocals of its elements and then transposing the resulting matrix..., *since multiplying something by a square diagonal matrix diag(v) is actually just calculating their* ***Hadamard*** *product... and here you want to get an identity matrix after multiplying it by its inverse*... ) D is not necessarily square tho; that's why it's 'pseudo'inverse.  
		- The pseudoinverse computes the solution that x* = A<sup>+</sup>y as it *projects* the solution x onto the *row space* of A, *which is orthogonal to the nullspace*!!!! And since it's literally ***orthogonal***, it's *geometrically the 'shortest' vector reaching the solution subspace* -> L2 norm is minimized!!!! This explains why: 
			-  When A has more columns than rows (more variables than independent equations; underdetermined), applying the pseudoinverse provides the solution x* = A<sup>+</sup>y with *minimal* L2/Euclidean norm ||x||<sup>2</sup> among all possible solutions of x.
		- Or, when A has more rows than columns (more equations than variable; overdetermined), it is possible for there to be no solution. In this case, using the pseudoinverse gives us the x for which *Ax is as close as possible to y* (the least squares solution) in terms of the minimized L2 norm ||Ax − y||<sup>2</sup>.
		- (I don't know whether there is a special meaning for substituting the notation of b (the known vector in the linear equation, as denoted earlier) to y in this section. Perhaps just to align with another set of terminology?)

### PCA 
Principal components analysis (PCA). (**It is essentially SVD applied to a mean-centered (standardized) data matrix.**) 

Re-read after reading Ch4 for numerical computation (minimized an expression by setting its gradient to zero).

The derivation afterwards was quite long)) but the key idea is here. The book also covers PCA in Ch5 on ML. 
![Screenshot 2025-04-17 at 09.15.51.png](https://img.picui.cn/free/2025/04/17/6800564c35565.png)



---

// I used to detest reading heavy math sections like this but atm after i've tried to keep notes and 'pause and ponder' parts that i couldn't understand... they don't seem that bad..!! and life feels hopeful...


## Chapter 3: Probability and Information Theory

// Skipped for some

- Probability theory is a mathematical framework for representing *uncertain* statements. In AI applications, used for: 
	1. The laws of probability tell us how AI systems should reason, so we design our algorithms to compute or approximate various expressions derived using probability theory. 
	2. Theoretically analyze the behavior of proposed AI systems.

### Probability Distributions
- Some notations:  ![1744852817460.png](https://img.picui.cn/free/2025/04/17/6800575666ed1.png)

#### Joint probability distribution and related concepts

- Marginal probability distribution: focusing on just a *subset* of the set of random variables with a given *joint* probability distribution. Computed by summing (for PMFs)/integrating (for PDFs); this summation gives its name 'marginal' (as it's usually carried out in the margin of a paper with the data chart).

- The **chain rule** of conditional probabilities:
	
	-  Used in *splitting a chunky joint probability distribution into many factors* (factorization; see below in **directed graphical models**). ![1744852843631.png](https://img.picui.cn/free/2025/04/17/6800576e3914d.png)

#### Common probability distributions
- **Bernoulli**: controlled by a single parameter φ ∈ [0, 1], which gives the probability of the random variable being equal to 1. -> P(x = 0) = 1 − φ
- **Multinoulli/categorical**: used to refer to a distribution over one single discrete variable but with k (finite) different states -> describing categories (eg. 'white', 'cat'). Since categories/labels are non-numeric, statistics like expectation or variance are irrelevant for multinoulli random variables. 
	- One notation to note here: 
	> The final, k-th state's probability is given by 1-1<sup>T</sup>p.
	
	I was confused about the notation '1<sup>T</sup>', but it turns out that this simply means a transposed *vector of ones*: [1, 1, ..., 1]). 
	p∈[0,1]<sup>k−1</sup>: A vector, its entries are probabilities for the first  k−1 categories. Also note that the possibility for the final/k-th state is no longer a free variable (and that's why there's the constraint 1<sup>T</sup> >= 0).

- **Gaussian/normal**: 
  
	- SUPER common for applications. In the absence of prior knowledge about what form a distribution over the real numbers should take, the normal distribution is a good default choice for two major reasons:
		1. Many distributions we wish to model are *truly* close to being normal distributions according to the **central limit theorem** (the *sum* of many independent random variables is approximately normally distributed. so in practice, many complicated systems can be modeled successfully as normally distributed noise, *even if the system can be decomposed into parts with more structured behavior.*)
		
		2. Out of all possible probability distributions *with the same variance (fixed)*, the normal distribution encodes the maximum amount of uncertainty (i.e., highest entropy). 
		    - Alongside calculating the entropy H<sub>Gaussian</sub> with the formula (which exceeds the entropy of any other distribution with the same variance), an intuition I find helpful is that the Gaussian spreads probability mass as 'evenly' around the mean as possible *given the variance constraint*. Other distributions (eg. Laplace below, sharper peaks and heavier tails; asymmetric) impose additional structure, which reduces uncertainty -> more predictable. 
		    - Because of this, we can think of the normal distribution as being the one that has to **assume the least amount of prior knowledge** into a model.
	- When the Gaussian distribution generalizes to R^n, it is known as the multivariate Gaussian distribution. Now the parameters μ is a *vector* and Σ gives the *covariance* matrix of the distribution.
		- *Isotropic Gaussian*: covariance matrix  Σ=σ<sup>2</sup>I, which means all pairwise covariances are zero -> all variables have exactly the same variance σ<sup>2</sup>. Requires only one parameter σ<sup>2</sup> instead of a full covariance matrix, so super simple and computationally efficient.
		- Percision matrix: the inverse of the covariance matrix, Λ = Σ<sup>-1</sup>. (In the context of a *uni*variate Gaussian distribution, precision is the reciprocal of the variance, τ = 1/σ². It quantifies how concentrated the distribution is around the mean).
		- Each multivariate Gaussian distribution (p(x | c = i), with *its own* mean µ^(i)^ and covariance matrix Σ^(i)^) can act as a component (c) and combine (via weighted sum, so *weight is another parameter for mixture models*) to form a Gaussian mixture model, elaborated below.
- **Exponential**: In DL, we often want to have a probability distribution with a sharp point at x = 0. And this distribution also uses the indicator function **1**<sub>x>=0</sub> mentioned in the 'Some notations' note section, where probability zero is assigned to all negative values of x. 
	- **Laplace**: Similarly, allows us to *place a sharp peak of probability mass* at an arbitrary point µ.
- **Dirac**: a bit abstract,,, often used to model a point mass or an infinitely sharp spike at a particular point so that the function is zero everywhere except at that particular point, and its integral over the entire real line is literally equal to one.
	- Used as *components* (again, see below for *mixture distributions*) in **empirical distributions**: each observed data point is assigned a probability mass of  1/n, represented by a delta function.
		- So why not just use a PMF for n discrete variables? -> Flexibility to use integration and other tools for continuous distributions. It also models discrete events within continuous time systems (eg. spikes in neuronal potentials!!!!)

- **Mixtures of distributions** :
	
	- **Gaussian mixture models**
	(GMMs) are considered as universal *approximators*, since by using a sufficient number of Gaussian components, we can approximate any smooth density function to any desired accuracy via the **weighted sum of these Gaussians** (similar to how NNs are universal approximators for functions...) The mixture coefficients, means, and covariances can be adjusted (through the weights assigned to each component) to fit the target density closely as the number of components increases.
	A GMM is a probabilistic model that assumes all the data points are modeled from *a mixture of several Gaussian distributions **with unknown parameters*** (these Gaussian distributions can describe *different clusters of data points*. 
	- The component identity variable c = i of a mixture model (in the former case, each multivariate Gaussian component) is a **latent variable** in the mixture model P(x) that we cannot observe directly (we can't juse 'see' P(c=i) in `P(x)`, can we?). They certainly are linked to x through the *joint distribution*, P(x) = P(x | c)P (c), but it is possible to describe P(x) without reference to the latent variable.
		- **Prior probability**: *α<sub>i</sub>* = P(c = i) given to each component i. “Prior” indicates that it expresses the model’s *beliefs about c before it has observed x*, representing the probability of that component being *active*. **This is just the weight of the i-th component for the mixture model! Since the mixture model is:**  ![1744852856642.png](https://img.picui.cn/free/2025/04/17/6800577b18918.png)
  
		- *Posterior probability*: P(c | x), since  it's *updating the belief* about that component *after observing x,* (the thingy calculated using Baye's Rule~~)![1744852879228.png](https://img.picui.cn/free/2025/04/17/6800579187cba.png)


#### Structured probabilitic models/graphical models
Machine learning algorithms often involve probability distributions over a *very large number* of random variables.
- Often, these probability distributions involve *direct* interactions between *relatively few variables*. So instead of using a single function to represent a probability distribution, we can *split* a probability distribution into many factors (for computational ease) *according to the **direct** interactions between them*.
- Eg. p(a, b, c) = p(a)p(b | a)p(c | b), where a influences the value of b and b influences the value of c, but that a and c are independent given b (so no direct interactions, as can be seen in the factorization).
- We can describe these *factorizations of probability distributions* using graphs; when we represent them with a graph, we call it a structured probabilistic model or graphical model. ૮₍˶ᵔ ᵕ ᵔ˶ ₎ა

> Throughout parts 1 and 2 of the book, structured probabilistic models will merely be used as a language to describe which **direct probabilistic relationships** different machine learning algorithms choose to represent.

**Directed** models use graphs with directed edges, and they represent factorizations into *conditional probability distributions*, as in the simple example above, and in this one: ![1744852900116.png](https://img.picui.cn/free/2025/04/17/680057a894525.png)

**Undirected** models use graphs with undirected edges, and they represent factorizations into a set of *functions*, but these functions are usually *not probability distributions of any kind; 'just functions'...* (denoted by φ<sup>(i)</sup>). Z is a normalizing constant to eventually obtain the normalized probability distribution p(x).![1744852918293.png](https://img.picui.cn/free/2025/04/17/680057b940cca.png)



### **Information theory** 
// i'm so extremely lazy here >:(

Information theory allows us to **quantify** the amount of uncertainty in a probability distribution.

The basic intuition behind information theory is that *learning that an unlikely event has occurred is more informative* than learning that a likely event has occurred. 
- So likely events should have low information content, and in the extreme case, events that are guaranteed to happen should have no information content whatsoever. 
- Independent events should have additive information.

#### Quantifying content of information 
The **self-information** of an event x = *x*: I(*x*) = − log P(*x*) 

> In the book, we always use log to mean the natural logarithm!!! 

Our definition of I(*x*) is therefore written in units of nats. One nat is the amount of information gained by observing **an event of probability 1/e** (since it's taking the *negative logarithm*).
- Other texts use base-2 logarithms and units called bits or shannons; information measured in bits is just a rescaling of information measured in nats.

**Shannon entropy** (*weighted average* (look at the expected value!!) of self-information):
![1744852951654.png](https://img.picui.cn/free/2025/04/17/6800583917e5b.png)
*When p = 0.5, the entropy is maximal, because the distribution is uniform over the two outcomes*.

#### Quantifying divergence of distributions

![1744852960157.png](https://img.picui.cn/free/2025/04/17/6800583a3ce14.png)

**Cross-entropy!!!!** (Comparing the formula below with the formula for Shannon entropy... it literally is *'cross'*-entropy.)
![1744852970219.png](https://img.picui.cn/free/2025/04/17/680058391183d.png) 



## Chapter 4: Numerical Computation
// as someone who has never learned multivariable calculus, this chapter helped so much more than i'd ever expect. before i've watched 3b1b's episode on gradient descent, which for me didn't provide a very convincing explanation (quite unusual for his videos...), and so was reading the topic on optimization from another much less in-depth DL book. so tbh it was such a surprise for this chapter to be so clear OwO!

---

**Numerical computation** typically refers to algorithms that solve mathematical problems by methods that **update** estimates of the solution via an **iterative process**, rather than analytically deriving a formula providing a symbolic expression for the correct solution (this is a spoiler for the following content, but an example can be finding x by directly solving the equation that *the gradient of this function is equal to 0* at the point x, which is possible in some cases).
- Common operations include: optimization (`argmin`/`argmax`) and solving systems of linear equations.

Real numbers can make performing continuous math difficult on computers, since they can't be represented *precisely* using a finite amount of memory.
- *Approximation error* when representing the number in a computer:
	- **Underflow**: numbers near zero are rounded to zero; causes undesirable consequences like ZeroDivisionError in subsequent operations. **Overflow**: numbers with large magnitude are approximated as ∞ or −∞; further arithmetic will usually change these infinite values into not-a-number values (placeholders).
		- One example of a function that must be stabilized (using specific tricks) against underflow and overflow is the *softmax* function. It's often used to predict the probabilities associated with a multinoulli (i.e., *categorical*) distribution.
	- Most ppl can simply rely on low level libraries that provide stable implementations without considering these numerical issues. Additionally, *Theano* is an example of a software package that automatically detects and stabilizes many common numerically unstable expressions that arise in the context of deep learning.

 **Conditioning**: how rapidly a function changes with respect to small changes in its inputs. Functions that change rapidly when their inputs are perturbed slightly can be problematic because *rounding errors* in the inputs can result in large changes in the output. 
- Condition number: for a square matrix that has eigendecomposition, it is the ratio of the magnitude of the largest and smallest eigenvalue (recall eigendecomposition, this reflects the largest divergence of change in the directions after the transformation is applied). When this number is large, a function that can be represented by the **inverse** of this matrix is particularly sensitive to numerical *rounding error* in the input. 
	- This sensitivity is an *intrinsic property of the matrix* itself, not the result of rounding error during matrix *inversion*. **Poorly conditioned** matrices *amplify* pre-existing errors when we multiply by the true matrix *inverse*.

### Gradient-Based Optimization
// This is where the cool part begins!! 

Most deep learning algorithms involve optimization of a function. The f(x) we want to minimize or maximize is called the objective function or criterion. 
- We usually phrase most optimization problems in terms of minimizing f(x). Maximization may be accomplished just by `argmin(−f(x))` (since we're trying to find the x* anyway).
	- When minimizing the objective function, we may also call it the cost function, loss function, or error function.

#### Why  'gradient-based' optimization?
- The loss functions we try to minimize usually have multiple inputs but only one single *scalar* output (so that different values of outputs are *comparable* and '*minimization*' would then make sense). MAE, MSE, cross-entropy, etc., are all scalar-valued to enable this comparability.
	- In this case (f : R^n -> R) where functions have multiple inputs, we must make use of **partial derivatives**, which measures how f(x) changes as *only* that particular variable x<sub>i</sub> increases with an infinitesimal step from the point x (calculating the partial derivative *with regard to* x<sub>i</sub> is not as complex as i thought,,, just regard the other terms that don't contain x<sub>i</sub> as constant terms and then differentiate as in single variable calculus). 
		- **The gradient generalizes the notion of derivative to the case where the derivative is *with respect to a vector*: the gradient of f is the *vector containing all of the partial derivatives* ⋆˙⟡. Element i of the gradient is the partial derivative of f with respect to x<sub>i</sub>.** In multiple dimensions, the critical points (here we want the function to move to the local minimum as the final outcome of minimization) are points where every element of the gradient is equal to zero.

- Since the derivative specifies how to scale a small *increase* in the input in order to obtain the corresponding change in the output, it's useful for minimizing a function since we can *reduce* f(x) by moving x in small steps *with the opposite sign of the derivative*, which already sounds pretty intuitive...!
	- A more *rigorous* demonstration of why **the direction in which f(x) decreases the fastest/steepest is just the opposite direction as the gradient (i.e., the negative gradient direction)** involves the concept of **directional derivative** *u*, a unit vector that quantifies the instantaneous rate of change of f(x) when moving in that direction: ![1744852981473.png](https://img.picui.cn/free/2025/04/17/6800583a4b670.png)

- A positive scalar determining the size of the step: the **learning rate** ϵ. So after a step taken x will now move by -ϵ∇<sub>x</sub>f(x) (the learning rate times the negative gradient of that point). 
- There are several ways to find ϵ. We don't want it to be too big in order to avoid overshooting the minimum (inadvertently going uphill in directions with *strong positive* curvature, think second derivative tests); whereas when the step size is too small, it won't make significant progress in other directions with less curvature- consider the second derivatives, which are collected together into the **Hessian matrix**. Because of this, gradient descent performs poorly when the Hessian matrix of our objective function has a large/poor condition number.

	>  Sometimes, we can solve for the step size ϵ that makes the directional derivative vanish. 
	
	- 'Making the directional derivative *vanish*': after taking the gradient descent step  x′=x−ϵ∇<sub>x</sub>f(x), the directional derivative of f at the new point x′ in the direction of the negative gradient becomes  *zero* (this ensures that moving steps further in the same direction no longer decreases the function, indicating a local minimum along that direction). Method for doing so is elaborated below, but the key idea is that we just adjust the value of ϵ *iteratively* until that directional derivative is zero, or solving ϵ* *analytically* for quadratic functions (...which are exactly the same as their 'approximations' when using *second-order* Taylor series expansion). This is known as the **exact line search**, which guarantees *optimal progress per step* but is computationally expensive (since iteration is most likely required).
	- **Line search**: just plug in several values of ϵ (heuristic?) to evaluate f(x') (i.e., f(x−ϵ∇<sub>x</sub>f(x)), and choose the ϵ that results in the smallest function value. 

- If the input *and output* are both vectors for our objective function (f : R^m -> R^n) rather than having a single *scalar* output as described above, it is conceivable that we'd have to collect all the partial derivatives into the **Jacobian *matrix*** instead of to the *gradient, which is sadly only a vector* that fails to house so many values. 
- Now the **Hessian matrix** (already mentioned earlier tho): this time back to the case where f : R^n -> R (otherwise even a Hessian *matrix* would prob fail to house so MANY values due to dimensionality..???). In this case the Hessian is the Jacobian of the *gradient*.
	- Due to symmetry (explained below), the eigenvectors (guaranteed to be orthogonal) specify **curvature** directions while the eigenvalues quantify curvature strength.
	- Then it's mostly related to *optimizing the step size* (which is literally the topic where we first mentioned the Hessian matrix; optimal step size depends on curvature strength and thus the *eigenvalues of the Hessian*):![1744853001694.png](https://img.picui.cn/free/2025/04/17/680058393077a.png)
![1744853011202.png](https://img.picui.cn/free/2025/04/17/6800583d784e9.png)

And: using additional information from the Hessian matrix (second derivatives; not just from the gradient/the first derivatives this time!!) to guide optimization:
![1744853019078.png](https://img.picui.cn/free/2025/04/17/6800583f2e77d.png)


#### Constrained optimization
While finding the minimized f(x), there might be some **constraints on the possible values of x**. Accessible points x (lie within the set S) are called feasible points.
- One simple approach to constrained optimization is simply to modify gradient descent taking the constraint into account.
- A more sophisticated approach is to design a different, *unconstrained* optimization problem whose solution can be converted into a solution to the original, constrained optimization problem ('just get rid of the constraints!!')
	- The Karush–Kuhn–Tucker (**KKT**) approach (generalizes the method of Lagrange multipliers λ, which only allows equality constraints but not **inequality constraints**. -> the KKT approach involves the *generalized* Lagrange function, which introduces another variable α, for inequality constraints, alongside λ for equality constraints; these two variables are called KKT multipliers) provides a very general solution to constrained optimization.
![1744853023567.png](https://img.picui.cn/free/2025/04/17/6800583f5b1de.png)

I watched a rlly cool [video](https://www.youtube.com/watch?v=GR4ff0dTLTw) that talks about how to understand the Lagrangian with only one equality constraint g(x) = 0, tho it's denoted as c(x) here.

Consider an objective function J(x) and an equality constraint c(x) = 0 (as denoted in the video screenshot below, in which we are looking at the intersection from a contour plot, and we are hiking yayyyy~~~). The super intuitive thing here is that **we are 'constrained' and thus can only walk on c(x) = 0, and we know that we have *reached the lowest point* of the objective function J(x) when our path is perpendicularly to its gradient ∇<sub>x</sub>J(x)**, since the negative gradient is just the *the direction of steepest decrease* (as demonstrated in the previous section "why  'gradient-based' optimization?"). So we want to find the location that satisfies this normal direction relationship.

- And... the **gradient of the constraint,  ∇<sub>x</sub>c(x),  just has the direction perpendicular to c(x) = 0!!!** (If you move  *along* the constraint set c(x)=0,  c(x) *doesn’t change!!!* (since c(x) = 0 everywhere on it -_-). Thus, the gradient of the set (i.e., the direction of steepest increase) must be perpendicular to the direction of the constraint set, which is the *direction of 'no change'*.
- And that's why this constrained optimization problem involving the equality constraint c(x) = 0 is converted to finding **the point x where the gradient of the objective function J(x) is parallel to the gradient of the constraint c(x)**!!!!! Since *gradients are just vectors (directional), we can scale one of them by a scalar* so that they would be equivalent (subtraction gives zero). *That's why we multiply the constraint by a constant, λ, the Lagrange multiplier*!!!!!!!! (just to avoid confusion, the author of the video says that the last row on the screenshot was written in that way (the term being added instead of subtracted) simply to make the equation look better since *λ is just a plain constant free to choose any sign, unlike the other KKT multiplier α* which must be >= 0 due to complementary slackness mentioned below)![Screenshot](https://img.picui.cn/free/2025/04/16/67ffb44315ea4.png)
That's why among the **KKT conditions**, the first one is that **the gradient of the generalized Lagrangian is zero** ('gives us the stationary points'!!). Then we just solve the equation involving the gradients by calculating the partial derivatives and stuff. Yayy!!
 
 Now also introducing the inequality constraints (and things get even more interesting...). Another KKT condition is the *complementary slackness* (α(h(x)) = 0): for every *inequality constraint*, either the constraint is _exactly_ enforced (α > 0 and *h(x) = 0*, the latter means that the inequality constraint is '**active**': the solution x after optimization (i.e., satisfies the first KKT condition that the gradient of the generalized Lagrangian is zero ->  the generalized Lagrangian function minimized) **lies exactly on the boundary** of the region of feasible points) or this particular inequality constraint just doesn’t matter at all  (*α = 0*, h(x) < 0) as we zero out the corresponding KKT multiplier α. 
 
<br>

 ![Screenshot 2025-04-16 at 23.36.44.png](https://img.picui.cn/free/2025/04/16/67ffce969ab20.png)

Now finally, concluding the KKT conditions (they're necessary conditions, but not always sufficient for a point to be optimal, **since in non-convex problems, the 'stationary point'** that satisfies KKT conditions can also be a local maximum/saddle point, so further verification is required). 
- The gradient of the generalized Lagrangian function (as shown in the screenshot right above) is zero.
-  All constraints on both x and the KKT multipliers are satisfied.
- The inequality constraints exhibit “complementary slackness”.

---
#### Applying constrained optimization: linear least squares
![1744857003660.png](https://img.picui.cn/free/2025/04/17/680067b1e5737.png)