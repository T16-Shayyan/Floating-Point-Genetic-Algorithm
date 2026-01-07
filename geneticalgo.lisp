;;; A SIMPLE GENETIC ALGORITHM OPERATING OVER FLOATING-POINT VECTORS


#|
1. Implement a very simple abstract high-level evolutionary computation framework

2. Implement a floating-point genetic algorithm

3. Test it on various objective functions and find good settings of parameters
which work well on those functions


|#





;;; Some utility Functions and Macros that might find to be useful

(defmacro while (test &rest body)
  "Repeatedly executes body as long as test returns true.  Then returns nil."
  `(loop while ,test do (progn ,@body)))

;;; Example usage
;;;
;;; (let ((x 0))
;;;    (while (< x 5)
;;;        (print x)
;;;        (incf x)))


(defun random? (&optional (prob 0.5))
  "Tosses a coin of prob probability of coming up heads,
then returns t if it's heads, else nil."
  (< (random 1.0) prob))

(defun generate-list (num function &optional no-duplicates)
  "Generates a list of size NUM, with each element created by
  (funcall FUNCTION).  If no-duplicates is t, then no duplicates
are permitted (FUNCTION is repeatedly called until a unique
new slot is created).  EQUALP is the default test used for duplicates."
  (let (bag)
    (while (< (length bag) num)
           (let ((candidate (funcall function)))
             (unless (and no-duplicates
                          (member candidate bag :test #'equalp))
               (push candidate bag))))
    bag))

;; hope this works right
(defun gaussian-random (mean variance)
  "Generates a random number under a gaussian distribution with the
given mean and variance (using the Box-Muller-Marsaglia method)"
  (let (x y (w 0))
    (while (not (and (< 0 w) (< w 1)))
           (setf x (- (random 2.0) 1.0))
           (setf y (- (random 2.0) 1.0))
           (setf w (+ (* x x) (* y y))))
    (+ mean (* x (sqrt variance) (sqrt (* -2 (/ (log w) w)))))))






;;;;;; TOP-LEVEL EVOLUTIONARY COMPUTATION FUNCTIONS


;;; TOURNAMENT SELECTION

;; is this a good setting?  Try tweaking it (any integer >= 2) and see
(defparameter *tournament-size* 4)
(defun tournament-select-one (population fitnesses)
  "Does one tournament selection and returns the selected individual."


  (let ((best-index (random (length population))))
    (dotimes (i (- *tournament-size* 1))
      (let ((next-index (random (length population))))
        (if (> (elt fitnesses next-index) (elt fitnesses best-index))
            (setf best-index next-index))))
    (elt population best-index)))



(defun tournament-selector (num population fitnesses)
  "Does NUM tournament selections, and puts them all in a list, then returns the list"


  (generate-list num (lambda () (tournament-select-one population fitnesses))))



;; I'm nice and am providing this for you.  :-)
(defun simple-printer (pop fitnesses)
  "Determines the individual in pop with the best (highest) fitness, then
prints that fitness and individual in a pleasing manner."
  (let (best-ind best-fit)
    (mapcar #'(lambda (ind fit)
                (when (or (not best-ind)
                          (< best-fit fit))
                  (setq best-ind ind)
                  (setq best-fit fit))) pop fitnesses)
    (format t "~%Best Individual of Generation...~%Fitness: ~a~%Individual:~a~%"
            best-fit best-ind)
    fitnesses))



(defun evolve (generations pop-size
               &key setup creator selector modifier evaluator printer)
  "Evolves for some number of GENERATIONS, creating a population of size
POP-SIZE, using various functions"

  
  ;; The functions passed in are as follows:
  ;;(SETUP)                     called at the beginning of evolution, to set up
  ;;                            global variables as necessary
  ;;(CREATOR)                   creates a random individual
  ;;(SELECTOR num pop fitneses) given a population and a list of corresponding fitnesses,
  ;;                            selects and returns NUM individuals as a list.
  ;;                            An individual may appear more than once in the list.
  ;;(MODIFIER ind1 ind2)        modifies individuals ind1 and ind2 by crossing them
  ;;                            over and mutating them.  Returns the two children
  ;;                            as a list: (child1 child2).  Nondestructive to
  ;;                            ind1 and ind2.
  ;;(PRINTER pop fitnesses)     prints the best individual in the population, plus
  ;;                            its fitness, and any other interesting statistics
  ;;                            you think interesting for that generation.
  ;;(EVALUATOR individual)      evaluates an individual, and returns its fitness.
  ;;Pop will be guaranteed to be a multiple of 2 in size.
  ;;
  ;; HIGHER FITNESSES ARE BETTER

  ;; your function should call PRINTER each generation, and also print out or the
  ;; best individual discovered over the whole run at the end, plus its fitness
  ;; and any other statistics you think might be nifty.

  
  (funcall setup)
  (let* ((population (generate-list pop-size creator))
         (best-individual nil))
    (dotimes (i generations)
      (format t "~%Generation ~a~%" (+ i 1))
      (funcall printer population (mapcar evaluator population))
      (mapcar (lambda (individual)
                (let ((fitness (funcall evaluator individual)))
                  (if (or (null best-individual) (> fitness (funcall evaluator best-individual)))
                      (setf best-individual individual))))
              population)
      (let ((new-gen nil))
        (dotimes (j (/ pop-size 2))
          (let* ((parent-a (first (funcall selector 1 population (mapcar evaluator population))))
                 (parent-b (first (funcall selector 1 population (mapcar evaluator population)))))
            (setf new-gen (append new-gen (funcall modifier parent-a parent-b)))))
        (setf population new-gen)))
    (format t "~%Best individual is ~a and its fitness is ~a"
            best-individual
            (funcall evaluator best-individual))))









;;;;;; FLOATING-POINT VECTOR GENETIC ALGORTITHM





(defparameter *float-vector-length* 20
  "The length of the vector individuals")
(defparameter *float-min* -5.12
  "The minimum legal value of a number in a vector")
(defparameter *float-max* 5.12
  "The maximum legal value of a number in a vector")

(defun float-vector-creator ()
  "Creates a floating-point-vector *float-vector-length* in size, filled with
UNIFORM random numbers in the range appropriate to the given problem"

  
  ;;;
  ;;; The numbers must be uniformly randomly chosen between *float-min* and
  ;;; *float-max*.  See the documentation for the RANDOM function.

  

  (generate-list *float-vector-length* (lambda ()
                                         (+ *float-min* (* (random 1.0)
                                                           (- *float-max* *float-min*))))))



;; I just made up these numbers, you'll probably need to tweak them
(defparameter *crossover-probability* 0.7
  "Per-gene probability of crossover in uniform crossover")
(defparameter *mutation-probability* 0.15
  "Per-gene probability of mutation in gaussian convolution")
(defparameter *mutation-variance* 0.05
  "Per-gene mutation variance in gaussian convolution")




;; to impement FLOAT-VECTOR-MODIFIER, the following two functions are
;; strongly reccommended.


(defun uniform-crossover (ind1 ind2)
"Performs uniform crossover on the two individuals, modifying them in place.
*crossover-probability* is the probability that any given allele will crossover. 
The individuals are guaranteed to be the same length.  Returns NIL."

  
  ;;; For crossover: use uniform crossover (Algorithm 25) in
  ;;;                Essentials of Metaheuristics

  
  (dotimes (i (length ind1))
    (if (random? *crossover-probability*) (rotatef (elt ind1 i) (elt ind2 i))))
  nil)


(defun gaussian-convolution (ind)
  "Performs gaussian convolution mutation on the individual, modifying it in place.
 Returns NIL."



  (dotimes (i (length ind))
    (if (random? *mutation-probability*)
        (let ((n (+ *float-max* 1.0)))  ;;;make sure n is out of bound to start with
          (while (or (< (+ (elt ind i) n) *float-min*) (> (+ (elt ind i) n) *float-max*))
                 (setf n (gaussian-random 0.0 *mutation-variance*)))
          (setf (elt ind i) (+ (elt ind i) n)))))
  nil)









(defun float-vector-modifier (ind1 ind2)
  "Copies and modifies ind1 and ind2 by crossing them over with a uniform crossover,
then mutates the children.  *crossover-probability* is the probability that any
given allele will crossover.  *mutation-probability* is the probability that any
given allele in a child will mutate.  Mutation does gaussian convolution on the allele."

    
    
    ;;; This function should first COPY the two individuals, then
    ;;; CROSS THEM OVER, then mutate the result using gaussian covolution,
    ;;; then return BOTH children together as a list (child1 child2)
    ;;;
    

  (let* ((child1 (copy-tree ind1))
         (child2 (copy-tree ind2)))
    (uniform-crossover child1 child2)
    (gaussian-convolution child1)
    (gaussian-convolution child2)
    (list child1 child2)))


;; you probably don't need to implement anything at all here
(defun float-vector-sum-setup ()
  "Does nothing.  Perhaps you might use this function to set
(ahem) various global variables which define the problem being evaluated
and the floating-point ranges involved, etc.  I dunno."
  )





;;; FITNESS EVALUATION FUNCTIONS

;;; I'm providing you with some classic objective functions.  See section 11.2.2 of
;;; Essentials of Metaheuristics for details on these functions.
;;;
;;; Many of these functions (sphere, rosenbrock, rastrigin, schwefel) are
;;; traditionally minimized rather than maximized.  We're assuming that higher
;;; values are "fitter" in this class, so I have taken the liberty of converting
;;; all the minimization functions into maximization functions by negating their
;;; outputs.  This means that you'll see a lot of negative values and that's fine;
;;; just remember that higher is always better.
;;;
;;; These functions also traditionally operate with different bounds on the
;;; minimum and maximum values of the numbers in the individuals' vectors.
;;; Let's assume that for all of these functions, these values can legally
;;; range from -5.12 to 5.12 inclusive.  One function (schwefel) normally goes from
;;; about -511 to +512, so if you look at the code you can see I'm multiplying
;;; the values by 100 to properly scale it so it now uses -5.12 to 5.12.


(defun sum-f (ind)
  "Performs the Sum objective function.  Assumes that ind is a list of floats"
  (reduce #'+ ind))

(defun step-f (ind)
  "Performs the Step objective function.  Assumes that ind is a list of floats"
  (+ (* 6 (length ind))
     (reduce #'+ (mapcar #'floor ind))))

(defun sphere-f (ind)
  "Performs the Sphere objective function.  Assumes that ind is a list of floats"
  (- (reduce #'+ (mapcar (lambda (x) (* x x)) ind))))

(defun rosenbrock-f (ind)
  "Performs the Rosenbrock objective function.  Assumes that ind is a list of floats"
  (- (reduce #'+ (mapcar (lambda (x x1)
                           (+ (* (- 1 x) (- 1 x))
                              (* 100 (- x1 (* x x)) (- x1 (* x x)))))
                         ind (rest ind)))))

(defun rastrigin-f (ind)
  "Performs the Rastrigin objective function.  Assumes that ind is a list of floats"
  (- (+ (* 10 (length ind))
        (reduce #'+ (mapcar (lambda (x) (- (* x x) (* 10 (cos (* 2 pi x)))))
                            ind)))))

(defun schwefel-f (ind)
  "Performs the Schwefel objective function.  Assumes that ind is a list of floats"
  (- (reduce #'+ (mapcar (lambda (x) (* (- x) (sin (sqrt (abs x)))))
                         (mapcar (lambda (x) (* x 100)) ind)))))




;;; an example way to fire up the GA.  If you've got it tuned right, it should quickly
;;; find individuals which are all very close to +5.12

#|
(evolve 50 1000
  :setup #'float-vector-sum-setup
  :creator #'float-vector-creator
  :selector #'tournament-selector
  :modifier #'float-vector-modifier
  :evaluator #'sum-f
  :printer #'simple-printer)
|#

#|
(evolve 50 1000
  :setup #'float-vector-sum-setup
  :creator #'float-vector-creator
  :selector #'tournament-selector
  :modifier #'float-vector-modifier
  :evaluator #'sphere-f
  :printer #'simple-printer)



(evolve 50 1000
  :setup #'float-vector-sum-setup
  :creator #'float-vector-creator
  :selector #'tournament-selector
  :modifier #'float-vector-modifier
  :evaluator #'rastrigin-f
  :printer #'simple-printer)
|#







;;;;;REPORT

;;;I started the project by reading through the given functions and compared what function 
;;;is going to map to what algorithm in the metaheuristics book. Then I started implementing
;;;the tournament-select-one function where I saw that the algorithm starts from 2 but then 
;;;realized that it's the same as iterating end - 1 times as we aren't using the index just 
;;;iterating that many times. Here I also learnt about elt function which helps find the ith 
;;;element in the list. After this I implemented tournament-selector function which was pretty 
;;;straight forward using the given generate-list function. After this I started the float-vector-creator 
;;;function which was algorithm 7 in the book. Here I learnt the unique way on how to use the 
;;;random function in lisp to find a random number between 2 bounds. After this was the 
;;;uniform-crossover function in which I learnt about the rotatef function which was interesting. 
;;;After that I did the gaussian-convolution function. Here I saw the while loop given as a macro. 
;;;In the book we used repeat until so while was the most appropriate use here. I was stuck 
;;;in how to declare n without initializing to go inside the while loop but found a work around 
;;;by declaring it as out of range. After this I did the float-vector-modifier where I learnt about 
;;;copy-tree to do a deep copy of a list from the cheat sheet. Finally I implemented the evolve function 
;;;which was algorithm 20 in the book. Here I learnt about how the passed in functions were as a pointer 
;;;so to run them I had to use funcall. I tried implementing it using the given professors functions but 
;;;got stuck inside so started using my own way in which I did dotimes instead of apply and walking step 
;;;by step with the algorithm making it easier as it was almost the same just with a fixed number of 
;;;generations. I also learned about the null function to check for nil. When experimenting around 
;;;with the parameters I noticed how when changing crossover probability for sum-f did not have much 
;;;difference but for a few others it slowed down the progress. So a middle probability would be the 
;;;best option in my opinion as it means it wont do too much crossover on the genes to produce noise 
;;;noe it would do less crossover to slow down improvements. When changing mutation probability a high 
;;;probability of around 0.5 caused too much fluctuation in the values that the overall fitness couldn't 
;;;reach the optimum. A really low mutation rate worked best for sum-f but wasn't the best option for 
;;;others like rastrigin-f. The best one I found was 0.08. Tournament size had the biggest effect out 
;;;of all. With a low tournament size of around 3 it was not reaching the best individual within the 
;;;50 generations and a high tournament value led to it reaching the best one really quickly which I 
;;;think is a classic trade off between exploitation and exploitation.















#|
Sum-f
Generation 1

Best Individual of Generation...
Fitness: 39.575306
Individual:(4.523773 4.9125805 3.722128 2.2933593 2.5855212 -1.6946485
            3.5684633 0.38613415 -0.47925186 1.7100916 4.329896 3.8634481
            -1.2414391 -1.6125062 1.6979518 1.5885401 -0.4097643 4.0817432
            4.92175 0.82753897)

Generation 2

Best Individual of Generation...
Fitness: 45.0343
Individual:(3.0026875 -1.7521143 1.5638685 2.3480506 4.6165724 4.4970884
            4.19724 -2.9663231 4.9213743 3.6143322 3.3652353 3.8054552 4.965211
            -1.7302489 4.5731373 2.62613 -0.21098042 1.3833556 2.8136744
            -0.5994444)

Generation 3

Best Individual of Generation...
Fitness: 55.948235
Individual:(0.33750725 4.0745945 3.6720266 3.9335861 1.7892823 2.6957812
            3.7700977 4.4943295 4.0702353 3.0644665 2.1643934 3.031145 4.402134
            4.000742 2.017849 4.9256687 -0.28567028 -0.17243528 0.96920156
            2.9933033)

Generation 4

Best Individual of Generation...
Fitness: 62.694412
Individual:(2.6323805 3.726919 3.541356 2.5453587 4.633749 2.903494 4.7567263
            3.5986204 -1.9071167 5.001405 3.4232168 4.2359877 4.730654
            0.91720724 3.6231723 5.036895 1.1932216 3.7887716 1.3190856
            2.9933033)

Generation 5

Best Individual of Generation...
Fitness: 66.48578
Individual:(4.0230923 1.9351015 4.273568 3.5859327 4.1892576 4.4835825
            1.2233152 4.0590463 2.8358226 3.3115606 4.689437 4.6497817
            4.9659357 2.6693187 1.4475427 1.5736756 5.0153437 3.1608448
            3.316187 1.0774302)

Generation 6

Best Individual of Generation...
Fitness: 77.201645
Individual:(4.69958 3.4935732 2.8350744 3.2246923 3.6305475 2.903494 2.057599
            3.985919 4.973031 4.108389 4.689437 4.528165 4.9659357 2.6693187
            5.052369 4.8010926 4.6129837 4.9950647 3.0664768 1.9089036)

Generation 7

Best Individual of Generation...
Fitness: 83.029396
Individual:(4.69958 3.4935732 4.6798573 2.5453587 3.5701122 4.5478115 4.7567263
            3.985919 4.973031 5.001405 4.689437 4.528165 4.9659357 2.6693187
            1.4475427 4.8010926 4.5534153 4.9950647 3.0664768 5.0595813)

Generation 8

Best Individual of Generation...
Fitness: 86.11496
Individual:(4.69958 3.4935732 4.775655 2.5453587 3.5701122 4.5478115 4.7567263
            3.985919 4.973031 5.061042 4.689437 4.528165 4.9659357 2.6886737
            4.3583155 4.8010926 4.5534153 4.9950647 3.0664768 5.0595813)

Generation 9

Best Individual of Generation...
Fitness: 88.88299
Individual:(4.69958 3.4935732 3.7831948 2.1315794 3.472867 4.5478115 4.7567263
            4.9242554 4.973031 5.001405 4.689437 4.528165 4.9659357 4.7297363
            4.154026 4.9302435 4.5534153 4.9950647 4.6496696 4.90328)

Generation 10

Best Individual of Generation...
Fitness: 89.83251
Individual:(4.69958 3.4935732 3.7831948 2.1315794 3.472867 4.5478115 4.7567263
            4.9242554 4.973031 5.001405 4.689437 4.528165 4.9659357 4.7297363
            5.1035366 4.9302435 4.5534153 4.9950647 4.6496696 4.90328)

Generation 11

Best Individual of Generation...
Fitness: 91.91296
Individual:(4.403347 4.7568555 4.1892033 3.9335861 5.0915356 4.4961514 4.862653
            4.8496885 4.973031 5.001405 4.353047 4.269001 4.9659357 2.6693187
            5.052369 4.8010926 4.6129837 4.9950647 4.6496696 4.9870243)

Generation 12

Best Individual of Generation...
Fitness: 92.4151
Individual:(4.314537 4.959749 4.2058077 5.0509295 4.633749 4.4961514 4.7567263
            4.9380093 4.10925 5.001405 4.689437 4.2359877 4.4965134 4.7297363
            5.052369 4.704906 4.9877367 4.9950647 4.6496696 3.4073637)

Generation 13

Best Individual of Generation...
Fitness: 93.695404
Individual:(4.403347 4.959749 4.2058077 3.9335861 5.0915356 4.804296 4.862653
            4.9380093 4.973031 5.001405 4.689437 4.528165 5.02372 3.9691553
            4.4311547 4.60943 4.8454585 5.075655 4.6496696 4.7001486)

Generation 14

Best Individual of Generation...
Fitness: 95.66494
Individual:(4.6419716 4.839632 4.775655 5.027407 5.0915356 4.584509 4.3739634
            4.9380093 4.973031 5.061042 4.744516 4.528165 4.9659357 3.9689999
            5.1035366 4.8022943 4.6129837 4.9950647 4.6496696 4.9870243)

Generation 15

Best Individual of Generation...
Fitness: 95.88671
Individual:(4.69958 4.9415154 4.9936876 3.9335861 5.0915356 4.4961514 4.7567263
            4.9242554 4.973031 4.9424825 4.539174 4.793108 4.9659357 4.7297363
            5.1035366 4.9133286 4.4575877 4.9950647 4.6496696 4.9870243)

Generation 16

Best Individual of Generation...
Fitness: 97.52614
Individual:(4.69958 4.959749 5.0055857 5.080808 5.010962 4.5478115 4.611887
            4.9380093 4.973031 5.001405 4.689437 4.4788632 4.9659357 4.7297363
            5.1035366 4.8010926 5.097553 4.9950647 4.776523 5.0595813)

Generation 17

Best Individual of Generation...
Fitness: 98.1671
Individual:(4.69958 4.959749 5.0055857 5.080808 5.010962 4.5478115 4.7567263
            4.9380093 4.973031 5.001405 4.689437 4.974983 4.9659357 4.7297363
            5.1035366 4.8010926 5.097553 4.9950647 4.776523 5.0595813)

Generation 18

Best Individual of Generation...
Fitness: 98.33746
Individual:(4.6419716 4.891547 5.0055857 5.027407 5.010962 4.4913626 5.076771
            4.9242554 4.973031 5.001405 4.689437 4.793108 4.9659357 5.0608706
            5.1035366 4.8010926 5.048028 4.9950647 4.776523 5.0595813)

Generation 19

Best Individual of Generation...
Fitness: 98.47048
Individual:(5.0742054 4.959749 5.0055857 5.027407 5.0915356 4.5478115 4.7567263
            4.9380093 4.973031 4.902983 4.689437 4.974983 4.9659357 4.7297363
            5.1035366 4.8010926 5.097553 4.9950647 4.776523 5.0595813)

Generation 20

Best Individual of Generation...
Fitness: 98.928635
Individual:(4.69958 4.871152 5.0055857 5.105762 5.010962 4.679856 5.0553703
            5.0539346 4.973031 5.001405 4.636864 4.974983 4.9659357 5.0608706
            5.1035366 4.8010926 5.097553 4.9950647 4.776523 5.0595813)

Generation 21

Best Individual of Generation...
Fitness: 99.140396
Individual:(4.7997875 4.891547 5.0055857 5.027407 5.010962 4.4913626 5.0591655
            4.9643054 4.973031 5.001405 4.7513247 4.933574 4.9659357 5.0608706
            5.1035366 5.101326 5.048028 4.9950647 4.8966064 5.0595813)

Generation 22

Best Individual of Generation...
Fitness: 99.47939
Individual:(5.0835266 4.959749 5.0676675 5.104118 5.0915356 4.5478115 5.076771
            4.9896727 4.973031 5.061042 4.744516 4.9130917 4.9659357 5.0608706
            5.016482 5.036895 4.9877367 4.9950647 4.776523 5.027356)

Generation 23

Best Individual of Generation...
Fitness: 99.68986
Individual:(5.0751743 4.959749 5.0055857 4.8512697 5.0915356 4.5478115 5.076771
            4.9896727 4.973031 5.061042 4.9550247 4.974983 4.9659357 5.0608706
            5.1035366 5.036895 5.097553 4.9950647 4.776523 5.0918317)

Generation 24

Best Individual of Generation...
Fitness: 99.90906
Individual:(5.0751743 4.959749 4.9583635 5.080808 5.0915356 4.5478115 5.076771
            4.9569445 4.973031 5.0919456 5.02946 4.974983 4.9659357 5.0608706
            5.1035366 4.9894543 5.1092567 4.9950647 4.776523 5.0918317)

Generation 25

Best Individual of Generation...
Fitness: 100.269775
Individual:(4.9545536 5.076503 5.0676675 5.1101875 5.0915356 5.0492506 5.076771
            5.0539346 4.973031 5.061042 4.7497263 4.9130917 4.9659357 5.0608706
            5.1035366 5.036895 5.097553 5.0238175 4.776523 5.027356)

Generation 26

Best Individual of Generation...
Fitness: 100.44184
Individual:(4.7880783 4.959749 5.0055857 5.0509295 5.112792 5.0492506 5.0553703
            5.0901675 4.9104886 5.061042 5.1163125 4.974983 4.890306 5.0608706
            5.1035366 5.0032625 5.097553 5.0730734 4.9789095 5.0595813)

Generation 27

Best Individual of Generation...
Fitness: 100.67408
Individual:(5.116323 4.994127 5.0933347 5.027407 5.0915356 4.8305163 5.076771
            5.0072017 4.973031 5.061042 5.106093 5.0943613 4.9659357 5.0608706
            5.1035366 5.036895 5.097553 5.0482373 4.776523 5.112784)

Generation 28

Best Individual of Generation...
Fitness: 100.97214
Individual:(5.116323 4.9415154 5.0676675 5.1101875 5.112792 4.908041 5.0553703
            5.0901675 4.973031 5.061042 5.1163125 4.974983 5.010992 5.0608706
            5.1035366 5.092409 5.097553 5.0730734 4.9789095 5.027356)

Generation 29

Best Individual of Generation...
Fitness: 101.058136
Individual:(5.116323 4.9415154 5.0676675 5.1101875 5.0915356 4.908041 5.0553703
            5.0901675 4.973031 5.061042 5.1163125 5.0822377 5.010992 5.0608706
            5.1035366 5.092409 5.097553 5.0730734 4.9789095 5.027356)

Generation 30

Best Individual of Generation...
Fitness: 101.19934
Individual:(5.116323 4.9415154 5.0676675 5.1101875 5.0915356 5.0492506
            5.0553703 5.0901675 4.973031 5.061042 5.1163125 5.0822377 5.010992
            5.0608706 5.1035366 5.092409 5.097553 5.0730734 4.9789095 5.027356)

Generation 31

Best Individual of Generation...
Fitness: 101.40149
Individual:(5.116323 4.9415154 5.0676675 5.1101875 5.0915356 5.0492506
            5.0553703 5.0901675 5.1171327 5.1092997 5.1163125 5.0822377
            5.0207815 5.0608706 5.1035366 5.092409 5.097553 5.0730734 4.9789095
            5.027356)

Generation 32

Best Individual of Generation...
Fitness: 101.51226
Individual:(5.116323 4.972586 5.0676675 5.1101875 5.0915356 5.0492506 5.0553703
            5.0901675 5.1171327 5.1092997 5.1163125 5.0822377 5.0207815
            5.0608706 5.1035366 5.092409 5.097553 5.0730734 5.058607 5.027356)

Generation 33

Best Individual of Generation...
Fitness: 101.71228
Individual:(5.116323 5.076503 5.0676675 5.10209 5.0915356 5.1011305 5.076771
            5.0901675 5.1171327 5.061042 5.1163125 5.116715 5.1179547 5.0608706
            5.1035366 5.019275 5.097553 5.0730734 5.0792685 5.027356)

Generation 34

Best Individual of Generation...
Fitness: 101.757484
Individual:(5.1050696 5.076503 5.0676675 5.105762 5.0915356 5.0492506 5.076771
            5.105156 5.1171327 5.1092997 5.1163125 5.0822377 5.113753 5.113169
            5.1035366 5.0368404 5.097553 5.0702662 5.060089 5.0595813)

Generation 35

Best Individual of Generation...
Fitness: 101.80716
Individual:(5.116323 5.076503 5.0676675 5.10209 5.0915356 5.1011305 5.076771
            5.105156 5.1171327 5.1092997 5.1163125 5.116715 5.1179547 5.113169
            5.1035366 5.019275 5.097553 5.0730734 5.058607 5.027356)

Generation 36

Best Individual of Generation...
Fitness: 101.833374
Individual:(5.116323 5.076503 5.0676675 5.10209 5.0915356 5.1011305 5.0553703
            5.0901675 5.1171327 5.1092997 5.1163125 5.116715 5.1179547
            5.0608706 5.1035366 5.092409 5.097553 5.0730734 5.100369 5.027356)

Generation 37

Best Individual of Generation...
Fitness: 101.90357
Individual:(5.116323 5.076503 5.0676675 5.10209 5.0915356 5.1011305 5.0663285
            5.0901675 5.1171327 5.1092997 5.1163125 5.116715 5.113753 5.0608706
            5.1035366 5.0929523 5.097553 5.0730734 5.0792685 5.111366)

Generation 38

Best Individual of Generation...
Fitness: 101.96899
Individual:(5.116323 5.0939555 5.051609 5.1101875 5.112792 5.0986094 5.110416
            5.105156 5.1171327 5.1092997 5.1163125 5.116715 5.1179547 5.113169
            5.1035366 5.019275 5.097553 5.0730734 5.0984993 5.0874257)

Generation 39

Best Individual of Generation...
Fitness: 101.96899
Individual:(5.116323 5.0939555 5.051609 5.1101875 5.112792 5.0986094 5.110416
            5.105156 5.1171327 5.1092997 5.1163125 5.116715 5.1179547 5.113169
            5.1035366 5.019275 5.097553 5.0730734 5.0984993 5.0874257)

Generation 40

Best Individual of Generation...
Fitness: 102.01909
Individual:(5.116323 5.0939555 5.0676675 5.10209 5.0915356 5.1011305 5.112519
            5.105156 5.1171327 5.1092997 5.1163125 5.116715 5.1179547 5.113169
            5.1035366 5.0925765 5.097553 5.0730734 5.058607 5.112784)

Generation 41

Best Individual of Generation...
Fitness: 102.07467
Individual:(5.116323 5.076503 5.1127944 5.1101875 5.0915356 5.1011305 5.076771
            5.0901675 5.1171327 5.1092997 5.1163125 5.116715 5.1179547 5.113169
            5.1035366 5.1039205 5.097553 5.0730734 5.1178074 5.112784)

Generation 42

Best Individual of Generation...
Fitness: 102.10861
Individual:(5.116323 5.115993 5.1127944 5.1101875 5.0915356 5.1011305 5.076771
            5.0901675 5.1171327 5.1092997 5.1163125 5.116715 5.1179547 5.113169
            5.097993 5.1039205 5.097553 5.0730734 5.1178074 5.112784)

Generation 43

Best Individual of Generation...
Fitness: 102.12789
Individual:(5.116323 5.0939555 5.1028543 5.1101875 5.112792 5.0986094 5.110416
            5.105156 5.1171327 5.1092997 5.1163125 5.116715 5.113753 5.113169
            5.1035366 5.107416 5.097553 5.0967975 5.0984993 5.0874257)

Generation 44

Best Individual of Generation...
Fitness: 102.14226
Individual:(5.116323 5.115993 5.1127944 5.1101875 5.0915356 5.1011305 5.110416
            5.0901675 5.1171327 5.1092997 5.1163125 5.116715 5.1179547 5.113169
            5.097993 5.1039205 5.097553 5.0730734 5.1178074 5.112784)

Generation 45

Best Individual of Generation...
Fitness: 102.17467
Individual:(5.116323 5.115993 5.1127944 5.1101875 5.0915356 5.1011305 5.110416
            5.0901675 5.1171327 5.1092997 5.1163125 5.116715 5.1179547 5.113169
            5.097993 5.1039205 5.097553 5.1054845 5.1178074 5.112784)

Generation 46

Best Individual of Generation...
Fitness: 102.16823
Individual:(5.116323 5.0939555 5.081566 5.1101875 5.112792 5.1011305 5.112519
            5.105156 5.1171327 5.1092997 5.1163125 5.116715 5.1179547 5.113169
            5.099762 5.1080866 5.097553 5.109444 5.1178074 5.111366)

Generation 47

Best Individual of Generation...
Fitness: 102.20142
Individual:(5.116323 5.115993 5.1028543 5.1175656 5.0915356 5.0986094 5.1100206
            5.0901675 5.1171327 5.1092997 5.1163125 5.116715 5.1179547 5.113169
            5.097993 5.1031966 5.1193795 5.116613 5.1178074 5.112784)

Generation 48

Best Individual of Generation...
Fitness: 102.20062
Individual:(5.116323 5.115993 5.1127944 5.1175656 5.112792 5.1011305 5.110847
            5.1036034 5.1171327 5.1092997 5.1163125 5.116715 5.1179547 5.113169
            5.1035366 5.092409 5.1193795 5.0730734 5.1178074 5.112784)

Generation 49

Best Individual of Generation...
Fitness: 102.22942
Individual:(5.116323 5.115993 5.1127944 5.1175656 5.112792 5.1011305 5.110416
            5.1036034 5.1164994 5.1092997 5.1163125 5.116715 5.1179547 5.113169
            5.1113105 5.1039205 5.097553 5.1054845 5.1178074 5.112784)

Generation 50

Best Individual of Generation...
Fitness: 102.24487
Individual:(5.116323 5.0939555 5.1127944 5.1047125 5.112792 5.119457 5.115058
            5.105156 5.1171327 5.1092997 5.1163125 5.116715 5.1179547 5.113169
            5.1035366 5.1039205 5.1193795 5.116613 5.1178074 5.112784)

Best individual is (5.116323 5.0939555 5.1127944 5.1047125 5.112792 5.119457
                    5.115058 5.105156 5.1171327 5.1092997 5.1163125 5.116715
                    5.1179547 5.113169 5.1035366 5.1039205 5.1193795 5.116613
                    5.1178074 5.112784) and its fitness is 102.24487






Sphere-f
Generation 1

Best Individual of Generation...
Fitness: -81.9872
Individual:(2.650497 -3.4585447 -1.6530406 0.61905766 1.3469496 -2.5473328
            1.9236803 2.4265637 2.7793508 -1.361692 1.2801623 -1.2948401
            1.267179 -1.2502015 -2.5920398 -0.108691216 0.77424574 0.9486499
            -0.6457362 -4.157258)

Generation 2

Best Individual of Generation...
Fitness: -59.991337
Individual:(-0.21918678 -1.6154406 -1.552119 1.8877964 2.6977353 -0.11145878
            1.2836232 -1.5292358 -1.3214526 -1.8683057 1.4921546 -2.2413561
            0.2426076 3.3548994 -0.2500453 -2.356826 -2.64911 -1.0014977
            1.2388344 1.0537963)

Generation 3

Best Individual of Generation...
Fitness: -43.426956
Individual:(-0.8854995 -1.3386328 2.0553112 0.033165026 0.936708 0.6945815
            0.8144655 0.21435289 2.2229764 1.2573571 -1.181328 -0.884603
            -0.7029028 -1.1410182 1.8822923 3.16928 -1.8164709 1.4420104
            -1.4990051 1.6871204)

Generation 4

Best Individual of Generation...
Fitness: -26.187717
Individual:(-1.3660657 1.4063282 0.028579235 2.42141 0.584857 -0.110669136
            -0.15392685 -0.1780591 -1.9215147 -0.27760887 0.5721655 0.030216217
            1.5114379 1.3276906 0.9945489 0.5583763 -0.5015893 -2.3728516
            -0.74117184 0.43935537)

Generation 5

Best Individual of Generation...
Fitness: -17.212294
Individual:(0.77434206 1.4545546 1.0100832 1.1107578 1.206903 -0.097614765
            0.0273633 0.6702881 -0.8599415 0.68416977 -0.5249877 0.59468776
            2.1318216 0.832407 -1.2889342 0.5583763 -0.3241992 0.474535
            -0.9195385 -0.32146144)

Generation 6

Best Individual of Generation...
Fitness: -11.613885
Individual:(-0.39394522 -0.27004957 -0.30481434 -0.7114601 0.6121792 0.8685436
            -0.6668105 0.9372987 -1.3214526 0.68416977 0.0013027191 0.013449162
            -0.25925273 -0.21540473 -0.70391226 -0.022126675 -0.1477519
            1.1975622 -1.8625464 -0.76499367)

Generation 7

Best Individual of Generation...
Fitness: -8.277292
Individual:(0.0015591085 -0.27004957 -0.13622196 -0.7114601 0.6121792 0.8685436
            -0.6668105 0.9372987 -0.6505941 0.68416977 0.0013027191 0.013449162
            -0.25925273 -0.21540473 -0.060135905 0.22521257 -0.1477519
            0.3028493 -1.8625464 -0.76499367)

Generation 8

Best Individual of Generation...
Fitness: -7.649887
Individual:(0.21914053 0.4912305 -0.45884362 0.007021427 0.92680645 0.10825347
            0.12698364 -0.42236176 -0.46001816 -0.40784848 -0.5249877
            -0.82659674 1.0039611 0.2543212 1.1149964 0.09176779 0.39935738
            0.9788451 -1.0549989 0.43935537)

Generation 9

Best Individual of Generation...
Fitness: -5.0153284
Individual:(0.0015591085 -0.23550558 -0.13622196 -0.7114601 0.4877426
            -0.11678844 -0.6668105 0.46168804 -0.25104108 0.68416977 -0.5249877
            0.013449162 -0.25925273 -0.38807964 -0.060135905 0.5583763
            -0.3241992 1.1975622 0.247406 -0.76499367)

Generation 10

Best Individual of Generation...
Fitness: -3.7725112
Individual:(0.0015591085 -0.5699088 -0.13622196 -0.6706028 0.4877426
            -0.11678844 -0.6668105 -0.8016138 0.092797406 0.68416977 0.35470724
            0.501338 -0.25925273 -0.35114384 -0.060135905 -0.15248135
            -0.41799656 -0.3748241 0.247406 0.43935537)

Generation 11

Best Individual of Generation...
Fitness: -2.7302787
Individual:(0.24264272 -0.28362894 -0.04948169 -0.256361 -0.764164 0.6725216
            0.12698364 0.18395996 0.1615853 -0.27760887 0.35470724 0.013449162
            0.7120948 -0.35114384 0.48768806 -0.013781906 -0.4263258
            -0.39634418 -0.0012097359 -0.012704805)

Generation 12

Best Individual of Generation...
Fitness: -1.1923586
Individual:(-0.3536811 -0.13009167 -0.13622196 -0.256361 0.31997824
            -0.097614765 0.0273633 -0.18560393 -0.032259554 -0.27760887
            0.0013027191 0.013449162 -0.14901972 0.1606884 0.48768806
            -0.09215927 -0.3241992 -0.39634418 0.21757656 0.36982536)

Generation 13

Best Individual of Generation...
Fitness: -1.0536872
Individual:(0.24264272 0.4912305 -0.04948169 -0.256361 0.31997824 -0.097614765
            0.12698364 0.16565122 0.1615853 0.24312115 0.0013027191 0.2477695
            0.4238357 -0.21540473 -0.060135905 -0.013781906 -0.29920667
            -0.17789252 -0.0012097359 0.17953253)

Generation 14

Best Individual of Generation...
Fitness: -0.85054564
Individual:(-0.18518433 0.08904329 0.18636455 0.007021427 0.09702435
            -0.097614765 0.02492614 0.21435289 0.1615853 -0.27760887
            -0.11416364 0.35314798 -0.25925273 0.058443546 -0.02205164
            0.09176779 -0.3241992 0.3028493 0.11455509 0.42145884)

Generation 15

Best Individual of Generation...
Fitness: -0.61815965
Individual:(0.21914053 0.08904329 0.18636455 0.007021427 0.31997824
            -0.097614765 -0.15392685 0.21435289 0.1615853 -0.27760887
            -0.11416364 0.10240644 -0.25925273 0.058443546 -0.02205164
            0.09176779 -0.3241992 0.09061527 0.11455509 -0.11489433)

Generation 16

Best Individual of Generation...
Fitness: -0.50296074
Individual:(0.008440018 -0.30708888 -0.13622196 -0.1464405 0.31997824
            -0.097614765 0.10608557 0.18395996 0.1615853 -0.27760887
            0.0013027191 0.013449162 -0.012992293 -0.13445675 -0.060135905
            -0.022126675 -0.22690661 0.04593754 -0.0057737525 0.17953253)

Generation 17

Best Individual of Generation...
Fitness: -0.36929128
Individual:(0.008440018 -0.008550122 0.028579235 -0.1464405 0.31997824
            -0.097614765 0.10608557 0.18395996 0.1615853 -0.27760887
            -0.09763104 0.013449162 -0.012992293 0.058443546 -0.060135905
            -0.022126675 -0.1477519 0.04593754 0.11455509 0.17953253)

Generation 18

Best Individual of Generation...
Fitness: -0.24649023
Individual:(0.008440018 -0.13009167 -0.048069246 0.007021427 0.14720787
            -0.097614765 0.12698364 -0.19032949 0.1615853 -0.04858996
            0.0013027191 0.030216217 -0.14901972 0.058443546 -0.044805035
            -0.022126675 0.0663278 0.04593754 0.21757656 0.17953253)

Generation 19

Best Individual of Generation...
Fitness: -0.1430285
Individual:(-0.07232542 -0.05532384 0.028579235 8.2304375e-4 0.11074829
            -0.02303996 0.10608557 -0.14287844 0.1615853 -0.04858996 0.11584753
            0.013449162 -0.012992293 -0.024143726 -0.060135905 -0.080495834
            0.046470135 0.04593754 -0.0057737525 0.17953253)

Generation 20

Best Individual of Generation...
Fitness: -0.1477497
Individual:(0.008440018 -0.13009167 -0.048069246 0.007021427 0.009087667
            -0.097614765 -0.09519369 -0.19032949 0.1615853 -0.06846161
            0.0013027191 0.10240644 -0.012992293 0.058443546 -0.060135905
            -0.022126675 -0.088018075 0.12875721 -0.0057737525 -0.012704805)

Generation 21

Best Individual of Generation...
Fitness: -0.103046685
Individual:(0.008440018 -0.13009167 -0.04948169 0.007021427 0.009087667
            -0.08119338 0.0273633 -0.0729095 6.31243e-4 0.062090516 0.15352742
            0.004064804 -0.07444771 0.064096145 0.037958864 -0.022126675
            -0.1477519 0.044386357 -0.014980301 0.08811995)

Generation 22

Best Individual of Generation...
Fitness: -0.04529326
Individual:(0.008440018 0.06319648 0.028579235 0.007021427 0.009087667
            -0.10685945 0.0273633 -0.08762024 6.31243e-4 0.062090516
            0.0013027191 0.004064804 -0.07444771 0.064096145 -0.045781884
            -0.022126675 0.046470135 0.044386357 -0.0057737525 -0.012704805)

Generation 23

Best Individual of Generation...
Fitness: -0.04836022
Individual:(0.008440018 0.06319648 0.028579235 0.007021427 0.009087667
            -0.097614765 0.0273633 -0.08762024 6.31243e-4 0.062090516
            0.0013027191 0.004064804 -0.07444771 0.058443546 -0.08801085
            -0.022126675 0.046470135 0.044386357 -0.0057737525 -0.012704805)

Generation 24

Best Individual of Generation...
Fitness: -0.03451687
Individual:(0.008440018 -0.008550122 0.028579235 0.007021427 0.0014900267
            -0.030337293 0.02492614 -0.0729095 -0.032259554 0.0012028664
            0.0013027191 0.013449162 -0.012992293 -0.024143726 -0.060135905
            -0.022126675 -0.09939394 0.04593754 -0.0012097359 0.09259354)

Generation 25

Best Individual of Generation...
Fitness: -0.028972387
Individual:(-0.0077626947 0.014459327 0.028579235 -0.030765083 0.0014900267
            -0.08119338 0.02492614 -0.01579322 6.31243e-4 0.062090516
            0.026408216 0.004064804 -0.07444771 -0.06428814 0.016595986
            -0.022126675 0.046470135 0.04593754 -0.0057737525 -0.012704805)

Generation 26

Best Individual of Generation...
Fitness: -0.016213922
Individual:(0.0015591085 -0.008550122 -0.04948169 0.007021427 0.01736784
            -0.030337293 0.003542509 0.04304202 6.31243e-4 0.062090516
            0.0013027191 -0.01234594 -0.012992293 -0.039098486 0.01064612
            -0.016383864 0.046470135 0.04593754 -0.0057737525 -0.012704805)

Generation 27

Best Individual of Generation...
Fitness: -0.014492955
Individual:(-0.0077626947 -0.008550122 -0.04948169 0.007021427 0.0014900267
            -0.003112793 0.0273633 -0.027745113 7.9483166e-4 0.0012028664
            0.0013027191 0.030216217 -0.012992293 -0.024143726 -0.060135905
            -0.013781906 0.01349923 0.06693012 -0.0057737525 -0.012704805)

Generation 28

Best Individual of Generation...
Fitness: -0.009867785
Individual:(0.008440018 0.014459327 0.028579235 8.2304375e-4 0.009087667
            -0.030337293 0.0273633 -0.01579322 6.31243e-4 0.0012028664
            0.0013027191 0.030216217 -0.012992293 -0.024143726 -0.045781884
            -0.022126675 0.046470135 -0.012682311 -0.0057737525 -0.012704805)

Generation 29

Best Individual of Generation...
Fitness: -0.0078043756
Individual:(0.008440018 -0.008550122 -2.782829e-4 0.007021427 0.0014900267
            0.018228061 0.0273633 -0.01579322 -0.018021412 0.0012028664
            0.0013027191 0.030216217 -0.012992293 -0.024143726 -0.02932678
            -0.03933263 0.01349923 0.012984864 0.03695921 -0.012704805)

Generation 30

Best Individual of Generation...
Fitness: -0.00573755
Individual:(0.008440018 -0.008550122 0.032669097 0.007021427 0.009087667
            -0.02303996 0.0273633 0.030354813 7.9483166e-4 0.0012028664
            0.0013027191 -0.027436972 -0.0010625981 -0.024143726 -0.0026100352
            -0.022126675 -6.6915154e-4 -0.012682311 -0.0057737525 -0.012704805)

Generation 31

Best Individual of Generation...
Fitness: -0.0036576074
Individual:(0.0015591085 -0.008550122 -0.015670683 0.007021427 0.009087667
            -0.003112793 0.02492614 -0.01579322 6.31243e-4 0.0012028664
            0.0013027191 0.013449162 -0.012992293 -0.032122634 0.002999913
            -0.022126675 -0.010441765 0.012984864 -0.0012097359 -0.012704805)

Generation 32

Best Individual of Generation...
Fitness: -0.0031000804
Individual:(0.0015591085 -0.008550122 -0.015670683 0.007021427 0.009087667
            -0.003112793 0.02492614 -0.01579322 6.31243e-4 0.0012028664
            0.0013027191 0.013449162 -0.012992293 -0.024143726 0.002999913
            -0.022126675 -6.6915154e-4 0.012984864 -0.0012097359 -0.012704805)

Generation 33

Best Individual of Generation...
Fitness: -0.003484488
Individual:(0.0015591085 -0.008550122 -0.009174075 0.007021427 0.0014900267
            -0.003112793 0.0273633 -0.03512472 6.31243e-4 0.0012028664
            0.0013027191 0.004064804 -0.012992293 -0.024143726 -0.0026100352
            0.012101628 -6.6915154e-4 -0.012682311 -0.0057737525 -0.012704805)

Generation 34

Best Individual of Generation...
Fitness: -0.0028194683
Individual:(0.008440018 0.014459327 0.028579235 8.2304375e-4 0.0014900267
            -0.003112793 0.003542509 0.0025943313 0.011340648 0.0012028664
            0.0013027191 0.013449162 -0.012992293 0.015567806 0.01064612
            -0.022126675 -6.6915154e-4 0.012984864 -0.0057737525 -0.012704805)

Generation 35

Best Individual of Generation...
Fitness: -0.0020607812
Individual:(0.008440018 -0.008550122 -0.022650614 8.2304375e-4 0.009087667
            -0.003112793 0.003542509 0.003115274 6.31243e-4 0.0012028664
            0.0013027191 -0.010723889 -0.012992293 -0.024143726 -0.0026100352
            0.0068780333 -6.6915154e-4 0.012984864 -0.0057737525 -0.012704805)

Generation 36

Best Individual of Generation...
Fitness: -0.0014234085
Individual:(0.008440018 -0.008550122 -0.008807085 8.2304375e-4 0.0014900267
            -0.012128338 0.003542509 0.0013805311 -0.018021412 0.0012028664
            0.0013027191 0.004064804 -0.012992293 -0.0057425275 -0.0026100352
            7.8142807e-4 0.01349923 -0.012682311 -0.0057737525 0.010344669)

Generation 37

Best Individual of Generation...
Fitness: -0.0010781803
Individual:(-0.009710078 -0.0066498453 -0.0077694952 8.2304375e-4 0.0014900267
            0.008877322 0.003542509 0.0013805311 6.31243e-4 0.0012028664
            0.0013027191 -0.01234594 -0.012992293 -0.0057425275 -0.0026100352
            7.8142807e-4 -0.010441765 0.012984864 -0.0057737525 0.010344669)

Generation 38

Best Individual of Generation...
Fitness: -9.3693286e-4
Individual:(0.0015591085 -0.008550122 0.0057629105 8.2304375e-4 0.009087667
            -0.003112793 0.003542509 0.0013805311 6.31243e-4 0.0012028664
            0.0013027191 -0.01234594 -0.012992293 -0.0018928312 0.01064612
            -0.0012495667 -6.6915154e-4 0.012984864 -0.0012097359 0.010344669)

Generation 39

Best Individual of Generation...
Fitness: -8.2728645e-4
Individual:(0.0015591085 -0.008550122 -0.0077694952 0.007021427 0.0014900267
            -0.003112793 0.003542509 0.0013805311 7.9483166e-4 0.0012028664
            0.0013027191 -0.01234594 -0.012992293 -0.0018928312 0.002999913
            7.8142807e-4 -6.6915154e-4 0.012984864 -0.0012097359 0.010344669)

Generation 40

Best Individual of Generation...
Fitness: -6.287272e-4
Individual:(0.0015591085 -0.008550122 -2.782829e-4 8.2304375e-4 0.0014900267
            -0.003112793 0.003542509 -0.01144661 6.31243e-4 0.0012028664
            0.0013027191 0.004064804 0.0015731901 -0.0018928312 -0.0026100352
            -4.6394393e-4 -6.6915154e-4 0.012984864 -0.0057737525 -0.012704805)

Generation 41

Best Individual of Generation...
Fitness: -6.6795223e-4
Individual:(0.0015591085 -0.008550122 -0.0077694952 8.2304375e-4 0.0014900267
            -0.003112793 0.003542509 -0.01144661 6.31243e-4 0.0012028664
            0.0013027191 0.004064804 0.0015731901 -0.0018928312 -0.0026100352
            0.004334869 -6.6915154e-4 -0.012682311 -0.0012097359 -0.012704805)

Generation 42

Best Individual of Generation...
Fitness: -4.732685e-4
Individual:(0.0015591085 -0.008550122 -2.782829e-4 8.2304375e-4 0.0014900267
            -0.003112793 0.003542509 0.003115274 6.31243e-4 0.0012028664
            0.0013027191 0.004064804 0.0015731901 9.598993e-4 -0.0026100352
            7.8142807e-4 -6.6915154e-4 0.012984864 -0.0012097359 -0.012704805)

Generation 43

Best Individual of Generation...
Fitness: -4.0494098e-4
Individual:(0.0015591085 0.0017892821 -2.782829e-4 8.2304375e-4 -2.0835106e-4
            -0.003112793 0.003542509 0.0067263544 7.9483166e-4 0.0012028664
            0.0013027191 0.004064804 0.010038685 -0.0018928312 0.002999913
            7.8142807e-4 -6.6915154e-4 -0.005782556 -0.0012097359 -0.012704805)

Generation 44

Best Individual of Generation...
Fitness: -3.841944e-4
Individual:(0.0015591085 -0.008550122 -2.782829e-4 0.007021427 0.0014900267
            -0.003112793 0.003542509 0.0025943313 6.31243e-4 0.0012028664
            0.0013027191 0.004064804 -0.0010625981 9.598993e-4 0.002999913
            -4.6394393e-4 -6.6915154e-4 -0.005782556 -0.0012097359 -0.012704805)

Generation 45

Best Individual of Generation...
Fitness: -3.24858e-4
Individual:(0.0015591085 0.0017892821 -2.782829e-4 8.2304375e-4 0.0014900267
            -0.003112793 0.003542509 0.0013805311 6.31243e-4 0.0056271553
            0.0013027191 0.004064804 0.0015731901 -0.0057425275 0.002999913
            7.8142807e-4 -6.6915154e-4 -0.005782556 -0.0012097359 -0.012704805)

Generation 46

Best Individual of Generation...
Fitness: -2.4600385e-4
Individual:(0.0015591085 0.0017892821 -2.782829e-4 8.2304375e-4 0.0014900267
            -0.003112793 0.003542509 0.0013805311 6.31243e-4 0.0012028664
            0.0013027191 0.004064804 -0.0023999335 -0.0018928312 0.002999913
            7.8142807e-4 -6.6915154e-4 -0.005782556 -0.0057737525 0.010344669)

Generation 47

Best Individual of Generation...
Fitness: -2.605546e-4
Individual:(0.0015591085 0.0020714058 -2.782829e-4 8.2304375e-4 0.0014900267
            -0.003112793 0.003542509 0.003115274 6.31243e-4 0.0012028664
            0.0013027191 -9.457208e-4 -0.0023999335 -0.0018928312 -0.0026100352
            -0.0012495667 -6.6915154e-4 -0.005782556 -0.0012097359 -0.012704805)

Generation 48

Best Individual of Generation...
Fitness: -1.7917885e-4
Individual:(0.0015591085 0.0020714058 -2.782829e-4 8.2304375e-4 0.0014900267
            -0.003112793 0.0033989195 0.0013805311 6.31243e-4 0.0012028664
            0.0013027191 0.004064804 -0.0010625981 -0.0018928312 0.002999913
            7.8142807e-4 -6.6915154e-4 0.0017389059 -0.0012097359 0.010344669)

Generation 49

Best Individual of Generation...
Fitness: -1.4103267e-4
Individual:(0.0015591085 0.0020714058 -2.782829e-4 8.2304375e-4 0.0014900267
            -0.003112793 0.003542509 0.003115274 7.9483166e-4 0.0012028664
            0.0013027191 0.004064804 -0.0023999335 -0.0018928312 0.002999913
            7.8142807e-4 -6.6915154e-4 -0.005782556 -0.0012097359 0.0049791764)

Generation 50

Best Individual of Generation...
Fitness: -1.3260492e-4
Individual:(0.0015591085 0.0020714058 -2.782829e-4 8.2304375e-4 0.0014900267
            -0.003112793 0.003542509 0.0013805311 6.31243e-4 0.0012028664
            0.0013027191 0.004064804 -0.0023999335 -0.0018928312 0.002999913
            -4.6394393e-4 -6.6915154e-4 -0.005782556 -0.0012097359 0.0049791764)

Best individual is (0.0015591085 0.0020714058 -2.782829e-4 8.2304375e-4
                    0.0014900267 -0.003112793 0.003542509 0.0013805311
                    6.31243e-4 0.0012028664 0.0013027191 0.004064804
                    -0.0023999335 -0.0018928312 0.002999913 -4.6394393e-4
                    -6.6915154e-4 -0.005782556 -0.0012097359 0.0049791764) and its fitness is -1.3260492e-4





Rastrigin-f
Generation 1

Best Individual of Generation...
Fitness: -228.36116823573d0
Individual:(-2.8826573 -4.0246606 -3.8513145 -4.15157 -4.02917 0.08504391
            -2.324968 -0.96247053 -1.963175 -0.47814465 0.1332922 -0.024599552
            1.1615806 -0.72747326 -0.08796406 -3.3811512 -0.042933464 2.9860716
            4.694659 -1.0193911)

Generation 2

Best Individual of Generation...
Fitness: -207.74255481165343d0
Individual:(-2.8826573 2.2183275 0.44817734 -4.15157 -4.02917 0.08504391
            2.799663 3.3940096 0.9699292 4.9253354 0.1332922 -0.024599552
            1.1615806 -0.72747326 -0.1310515 2.7388487 -0.042933464 2.9860716
            0.0735569 -1.0193911)

Generation 3

Best Individual of Generation...
Fitness: -197.91863524296386d0
Individual:(0.4520898 0.89652824 0.97887707 -0.24844742 2.9693384 -0.1791172
            1.030541 -2.910907 -0.15306044 1.7998873 0.9472289 -3.7011669
            -1.753026 -1.8157713 2.7660413 3.0007448 3.1269379 -3.9594016
            2.8519425 -1.1769457)

Generation 4

Best Individual of Generation...
Fitness: -160.54159098784132d0
Individual:(2.846395 0.045346737 0.7328944 1.2512628 0.11260925 1.0480089
            -4.890776 -0.178979 -2.2020483 -2.9145997 0.92840195 1.8619838
            1.1393055 -0.93495464 -0.85488415 -0.122323036 4.826125 0.91090596
            0.03544569 2.8410106)

Generation 5

Best Individual of Generation...
Fitness: -124.41605717818379d0
Individual:(-2.2167945 1.0181025 2.0805006 -0.15661949 -0.80711454 4.0582657
            0.94180536 2.9209092 1.9938651 1.1989778 -0.91529 1.9616259
            -0.15167122 -0.18627071 1.8822923 -0.17058945 0.20737076
            -0.09224552 2.2354035 -1.0360504)

Generation 6

Best Individual of Generation...
Fitness: -114.69369537207959d0
Individual:(1.1220739 2.116979 -1.0889297 0.82871926 -0.9691882 2.9779549
            0.03130298 -0.88233453 -2.8569481 2.1782017 -0.87058055 -0.9376378
            0.26973268 -0.0025683641 3.9699523 -0.01323159 0.80858755 0.0696125
            1.2416258 -1.1795189)

Generation 7

Best Individual of Generation...
Fitness: -94.88495241955934d0
Individual:(1.9317737 -0.0066813454 1.1431007 1.0966369 0.21139634 2.0231128
            3.131253 -2.1201427 1.1862012 2.0226657 0.01701641 0.038595583
            1.1002016 -0.043241464 0.8245547 1.0347056 0.80462533 -0.081618786
            1.2896719 0.85983276)

Generation 8

Best Individual of Generation...
Fitness: -94.51645412720667d0
Individual:(1.9317737 -0.2641333 1.1431007 1.2621655 0.21139634 3.1033227
            -0.02409711 1.0038099 1.0784143 1.9480644 0.01701641 0.038595583
            1.1002016 -0.043241464 0.8245547 1.0347056 -0.6484051 -0.081618786
            2.0350523 0.9995551)

Generation 9

Best Individual of Generation...
Fitness: -87.75390841183584d0
Individual:(1.9317737 -0.017177626 0.82906723 1.9264771 1.87139 1.9633354
            -0.08649653 -1.0937304 -2.8569481 -2.0261 0.1517916 0.038595583
            1.1002016 0.9997387 0.8245547 0.7128053 0.80462533 -0.081618786
            1.0669456 0.12464407)

Generation 10

Best Individual of Generation...
Fitness: -64.70703602142379d0
Individual:(1.9317737 -0.9201478 0.028579235 0.014725208 -0.046016008
            -0.021189898 -0.02409711 1.0038099 1.0847692 1.9480644 2.116836
            1.8619838 1.1002016 -0.06722364 -0.79098284 1.0347056 -0.6484051
            0.017629884 2.0350523 -0.9933543)

Generation 11

Best Individual of Generation...
Fitness: -62.76412087819813d0
Individual:(2.846395 0.122720346 0.028579235 0.07494892 -0.12375367 2.9260292
            -0.10636174 1.1376989 0.0568285 -1.9926379 0.1332922 1.9412838
            1.0125452 -0.93495464 0.9670371 -0.01323159 0.015434742
            -0.0074120983 1.042353 2.988545)

Generation 12

Best Individual of Generation...
Fitness: -53.30021923228327d0
Individual:(0.04825872 1.0719337 -0.11284854 0.9289245 0.027044736 0.013785146
            -0.10636174 0.9247784 2.1777754 1.0308076 0.1332922 0.01658018
            1.021077 -1.9922705 0.19727154 -0.09473239 0.09640787 -2.1030614
            1.042353 -0.87167263)

Generation 13

Best Individual of Generation...
Fitness: -55.71667728136492d0
Individual:(-2.0214734 1.0181025 0.07092621 1.0310044 -0.08321427 1.1002846
            -0.0018677711 0.87962574 1.8745226 0.96148694 -0.9050645
            0.058598027 2.8538842 -0.93844503 -0.10666911 1.1416094 -1.043047
            2.0168302 1.0482159 -1.0360504)

Generation 14

Best Individual of Generation...
Fitness: -49.43726329576154d0
Individual:(0.90730083 0.009820372 -1.0450888 0.8563472 0.020053072 1.0480089
            0.94180536 1.1133368 -0.9556956 -0.9939458 2.2142637 1.9841101
            1.1002016 -0.11087486 0.048787266 0.19123581 -0.90682745
            -0.034906194 -0.005263239 -0.00840725)

Generation 15

Best Individual of Generation...
Fitness: -40.872472867926916d0
Individual:(1.9317737 -0.017177626 2.0287538 0.07494892 -0.10297316 0.9485949
            1.0038544 1.0124226 -1.963175 -0.1529728 0.1332922 1.0617778
            1.0150204 -0.93495464 -0.058899775 0.9706542 0.015434742
            -0.0074120983 1.0378913 0.83334345)

Generation 16

Best Individual of Generation...
Fitness: -45.58057332413591d0
Individual:(0.04662609 2.0840287 -0.009231091 -1.0578723 -0.9691882 0.12701347
            0.9606267 -1.0272406 -1.2195172 -0.007786751 -0.9050645 0.12600179
            1.0256317 0.0619393 2.1190336 -0.067846775 -0.10548303 0.9455263
            -1.0047168 -1.0360504)

Generation 17

Best Individual of Generation...
Fitness: -35.28310985427615d0
Individual:(0.008440018 1.0152601 0.998422 1.974325 2.027046 1.0480089
            0.92643464 0.8910766 0.14275318 -2.0261 0.07679233 -0.024599552
            -0.9519914 -0.9488023 0.9716957 -0.043711245 1.0222584 -0.034906194
            -0.10132207 0.08871446)

Generation 18

Best Individual of Generation...
Fitness: -37.46874358564642d0
Individual:(0.07183838 1.0778596 -0.009231091 0.014725208 -0.9772337 0.8590412
            0.077762544 1.1159011 1.0763829 0.93124276 1.0800337 -0.9376378
            0.96995133 0.097911 1.8748511 0.007642174 -1.0773535 0.96856236
            0.08429017 1.0325437)

Generation 19

Best Individual of Generation...
Fitness: -32.52248353444827d0
Individual:(-0.024940372 2.033364 0.81487036 0.01701413 -0.94035983 1.0879662
            -0.027561933 1.0038099 0.08896144 0.03835881 -0.07196111
            0.080804035 -0.9519914 -0.93495464 0.9716957 -1.0372944 -0.07770407
            -0.041166637 -1.0174738 -1.1132866)

Generation 20

Best Individual of Generation...
Fitness: -22.706617415862496d0
Individual:(-0.024940372 1.9800415 -0.9830947 0.96457005 -0.031239433
            0.0042093173 0.07743963 0.03326539 0.08896144 0.98918486 -0.91529
            -0.024710134 -0.07721762 1.0017773 0.9716957 -0.043711245
            -0.90682745 -0.09224552 1.010663 -1.0193911)

Generation 21

Best Individual of Generation...
Fitness: -29.355706142608824d0
Individual:(-0.024940372 -0.017177626 0.998422 0.01701413 -0.89373577 1.0879662
            -0.10636174 1.0038099 0.08896144 0.03835881 -0.07196111 0.058598027
            1.1044377 -0.93495464 -0.042780876 -1.0372944 1.9753556
            -0.041166637 -1.0174738 -1.1132866)

Generation 22

Best Individual of Generation...
Fitness: -29.40850659692623d0
Individual:(-0.024940372 -0.017177626 -0.011186957 0.9412938 -1.1257938
            0.9072216 -0.10636174 1.0038099 0.08896144 1.0005589 0.99313766
            0.12303504 -0.047130466 1.0190477 -0.042780876 -0.97785187
            -0.90682745 1.0245423 -1.0174738 -1.1132866)

Generation 23

Best Individual of Generation...
Fitness: -29.57518980877674d0
Individual:(0.045689836 0.13187994 1.0264715 2.0391564 0.05728157 1.0357614
            0.07743963 0.98769695 0.08896144 0.9414134 1.007375 1.0435939
            1.1307851 -0.9488023 -0.009413585 0.01720506 -0.94397724
            0.052588075 0.062322643 -0.93176)

Generation 24

Best Individual of Generation...
Fitness: -29.710214778203408d0
Individual:(1.0430627 -0.017177626 0.22026765 0.9289245 -0.06481039 1.0357614
            -0.045588154 0.93553925 -0.013099007 0.03835881 0.045217957
            -0.024710134 1.0150204 -0.0025683641 -0.10666911 0.015298785
            0.951064 1.1055284 1.9685321 0.996309)

Generation 25

Best Individual of Generation...
Fitness: -24.286149164124765d0
Individual:(0.96328855 0.009820372 0.96463287 0.07379127 -0.06481039 0.91274744
            -0.043882005 0.9247784 -0.032322884 -0.009026825 -0.09006792
            0.06789678 -0.042802453 -0.9779414 -0.16785271 -1.0372944
            0.10423644 -0.041166637 -1.0174738 -1.0360504)

Generation 26

Best Individual of Generation...
Fitness: -22.38649144521858d0
Individual:(-0.00732287 0.91347826 0.9985118 0.9595113 -0.92114276 -0.051916778
            -0.056321546 -0.08742466 0.08896144 1.0926738 -0.052191615
            -0.024710134 0.9887099 -0.0025683641 0.9619235 0.07451671
            -1.0773535 0.017629884 1.010663 0.9479404)

Generation 27

Best Individual of Generation...
Fitness: -27.05812004357503d0
Individual:(-0.024940372 0.009820372 -0.11649616 0.01701413 0.14199677
            -0.1735731 0.9606267 1.0199238 -1.0632367 -0.055799037 -0.04476737
            -0.024710134 -0.9519914 0.011063695 0.031342633 0.021467265
            0.8488784 1.0406163 0.025312766 1.052468)

Generation 28

Best Individual of Generation...
Fitness: -26.390918989992997d0
Individual:(0.08284947 -0.1109056 0.1271527 0.07494892 0.07799578 -0.013284557
            0.8966481 0.97131234 0.08896144 0.94440234 -0.10163863 0.06455566
            0.9026206 -1.0547241 0.031342633 -0.015598752 0.04129698 0.05077967
            1.0260801 -0.9933543)

Generation 29

Best Individual of Generation...
Fitness: -25.726286261747845d0
Individual:(-0.014680788 -0.05221384 0.96463287 -1.0854332 0.024464704
            0.9562036 0.11051842 -0.009975687 0.95028055 1.0517144 -0.9526571
            0.044780977 1.025917 -0.053374708 0.04790011 0.1152178 -1.043047
            1.0245423 0.98753643 0.85655075)

Generation 30

Best Individual of Generation...
Fitness: -23.4575018780439d0
Individual:(-0.014680788 -0.05221384 0.9688343 -1.0854332 0.024464704 0.9562036
            0.11051842 -0.009975687 0.95028055 1.0796155 -0.9526571 0.044780977
            1.0681313 -0.053374708 0.02679596 0.1152178 -1.043047 1.0245423
            0.98753643 -0.9933543)

Generation 31

Best Individual of Generation...
Fitness: -21.107302382460517d0
Individual:(1.003444 0.9978669 0.04668607 0.01701413 0.024464704 1.0699137
            1.0288811 -0.0077789277 -0.054208845 0.04464049 -0.9526571
            0.044780977 -0.09619354 -0.0025683641 0.02679596 0.1152178
            1.0391169 1.0245423 0.85786676 -0.9933543)

Generation 32

Best Individual of Generation...
Fitness: -24.87822709889042d0
Individual:(0.040858343 0.019924022 0.8230251 -0.013239333 0.057950646
            1.0453991 0.020423107 -6.2985346e-4 1.0797327 1.0481484
            -0.037226006 -1.0298142 1.0390711 -1.0477176 0.9493637 0.13123916
            -1.043047 0.023166128 0.015892442 -1.0360504)

Generation 33

Best Individual of Generation...
Fitness: -23.0993560035638d0
Individual:(0.04299198 -0.05221384 0.9688343 0.035623066 0.13960849 -0.06904389
            0.07743963 1.0678147 0.97127765 -0.06080128 0.02842135 -0.08381194
            1.0390711 0.0619393 0.027127132 -0.00920704 -1.043047 1.0245423
            0.8804607 -1.0360504)

Generation 34

Best Individual of Generation...
Fitness: -11.066855746494241d0
Individual:(-0.00634331 1.0703069 -0.05188456 0.047636747 -0.06481039 1.027252
            -0.0018677711 0.010922372 0.0013831705 0.008321241 -0.059283882
            0.058598027 1.0150204 0.05417285 0.0022326708 -0.037636265
            1.0391169 -0.03885278 -0.039448544 0.056865208)

Generation 35

Best Individual of Generation...
Fitness: -17.507422823810543d0
Individual:(-0.00634331 -0.061841164 1.136497 0.047636747 0.1568115 1.027252
            -0.0018677711 0.010922372 -0.048355028 0.008321241 -0.059283882
            0.058598027 1.0150204 -0.009312995 0.0022326708 -0.037636265
            1.0391169 -0.03885278 -0.039448544 0.056865208)

Generation 36

Best Individual of Generation...
Fitness: -19.335097680893654d0
Individual:(0.08973096 -0.017177626 1.0119203 0.01701413 -0.06481039 0.9562036
            -0.017916873 -0.012347853 -0.12013312 0.061844815 -0.019041374
            -0.02995839 1.021077 -0.1397783 0.031342633 0.94788057 -0.9564562
            1.0269998 0.06274521 0.9814753)

Generation 37

Best Individual of Generation...
Fitness: -16.70693374489622d0
Individual:(-0.031174213 -0.022387289 0.062373064 0.027343815 -0.015862562
            0.9465201 0.07743963 1.0038099 0.10563946 0.03835881 -0.026856361
            -0.02995839 0.99749297 0.9301616 0.031342633 -0.015598752
            -0.9564562 0.016038395 0.119560584 0.9235937)

Generation 38

Best Individual of Generation...
Fitness: -19.921481570700877d0
Individual:(0.009726286 1.0008961 1.0736897 1.0143522 0.057950646 1.1016884
            0.027686082 -0.009975687 -0.054208845 0.9414134 0.02842135
            0.0073002726 1.0459858 -1.0477176 -1.0121008 -0.07593779 -1.0773535
            -0.0528862 -0.07877028 0.003027141)

Generation 39

Best Individual of Generation...
Fitness: -17.319101294910865d0
Individual:(0.025899455 -0.030672729 1.0647653 -0.046590805 1.0130925
            0.87049055 -0.01805681 0.03326539 0.97127765 -0.014814444 1.0544753
            -0.027900003 1.021077 0.060601786 1.0385971 -0.046623573
            -0.070861936 0.010570623 -0.07877028 0.052199215)

Generation 40

Best Individual of Generation...
Fitness: -16.565733553428174d0
Individual:(0.035890102 0.009820372 0.94516295 -0.017203156 0.05728157
            1.0566928 0.07743963 1.0038099 0.0013831705 0.03835881 -0.037226006
            -0.025616959 1.0837201 0.042284578 0.027127132 -0.081986524
            1.0058391 0.97682446 0.09683606 1.0329785)

Generation 41

Best Individual of Generation...
Fitness: -16.459586001831184d0
Individual:(0.035890102 0.009820372 0.9396778 0.01892431 0.05728157 1.0566928
            0.07743963 1.0678147 0.0013831705 0.03835881 -0.037226006
            -0.025616959 1.0837201 0.042284578 0.027127132 -0.081986524
            1.0058391 0.97682446 0.09683606 0.011670751)

Generation 42

Best Individual of Generation...
Fitness: -15.620323935553813d0
Individual:(0.08658476 0.9833345 -0.026973259 0.0100725815 -0.06669189
            1.0104892 -0.045588154 0.020227194 0.087152794 0.008321241
            -0.052191615 0.01658018 0.96902674 -0.0025683641 0.08661603
            -0.050787844 -1.0221573 0.97682446 0.940473 1.052214)

Generation 43

Best Individual of Generation...
Fitness: -12.790682782167153d0
Individual:(0.0084840655 -0.025630433 0.998422 0.040843338 0.100824475 1.05119
            -0.041448653 1.0678147 -1.0090002 -0.014814444 -0.0026038624
            0.04103028 1.0015874 0.005277492 0.031342633 -0.03269638
            0.037161008 0.97682446 -0.033675075 0.99653095)

Generation 44

Best Individual of Generation...
Fitness: -14.444924959646102d0
Individual:(0.047211036 -0.017177626 0.9396778 0.0100725815 0.14504576
            0.9956922 1.0288811 -0.0077789277 -0.0052362457 0.008321241
            0.9078749 -0.02995839 1.0150204 -0.0025683641 -0.052352488
            0.06289636 0.0631492 0.017629884 0.03357585 -0.023089569)

Generation 45

Best Individual of Generation...
Fitness: -15.116297020030629d0
Individual:(0.047211036 0.009820372 0.9396778 0.0100725815 0.14504576 0.9956922
            0.12256484 -0.0077789277 -0.048471708 0.008321241 -0.061688855
            -0.02995839 1.021077 -0.0025683641 -0.052352488 0.06289636
            -0.97510326 0.017629884 0.03357585 -0.023089569)

Generation 46

Best Individual of Generation...
Fitness: -20.572006049603146d0
Individual:(0.025899455 -0.12869588 0.06875122 -0.017203156 0.025553819
            0.956589 -0.04511685 -0.11660926 0.0013831705 -0.040489443
            0.09991343 -0.041695282 0.92301697 0.0561565 0.06975038 0.10103138
            0.07205282 0.017629884 0.8999032 0.014776047)

Generation 47

Best Individual of Generation...
Fitness: -14.946497694842037d0
Individual:(1.037641 -0.061841164 -0.05188456 0.040843338 -0.10554719 1.0198603
            0.053067602 -0.09746022 0.07239838 -0.02103658 -0.029510247
            -0.0058644563 -1.0932993 0.052569844 0.027127132 0.01720506
            0.04471883 0.032696776 0.042604182 0.049898244)

Generation 48

Best Individual of Generation...
Fitness: -17.22578063891146d0
Individual:(0.025899455 -0.061841164 -0.024642192 -0.017203156 -0.008052326
            0.99176526 1.032023 -0.11660926 0.10915646 0.082566865 0.025043927
            0.066990495 0.9801117 -1.0477176 0.019534871 0.01720506 1.0222584
            0.010570623 1.025208 0.09258809)

Generation 49

Best Individual of Generation...
Fitness: -16.53441933233549d0
Individual:(-0.043030433 -0.075563565 0.04133673 1.007473 0.00612285 0.99818027
            0.06538445 1.0678147 -0.049351893 1.1099911 0.06483887 0.03221041
            -0.03094349 0.052569844 0.027127132 -0.029614702 -0.95614463
            0.9931627 -0.009472609 0.9479404)

Generation 50

Best Individual of Generation...
Fitness: -16.086013687832946d0
Individual:(-0.07567478 0.12092696 -0.052551456 -0.046590805 0.1271266
            0.9929385 0.06218977 -0.009975687 -0.0688947 -0.034313448
            0.07045816 0.06628113 1.0219014 0.012846938 0.019534871 0.030896591
            -0.004870534 1.0346125 -0.038392525 0.044860546)

Best individual is (-0.00634331 1.0703069 -0.05188456 0.047636747 -0.06481039
                    1.027252 -0.0018677711 0.010922372 0.0013831705 0.008321241
                    -0.059283882 0.058598027 1.0150204 0.05417285 0.0022326708
                    -0.037636265 1.0391169 -0.03885278 -0.039448544 0.056865208) and its fitness is -11.066855746494241d0


|#


