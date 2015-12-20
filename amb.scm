;;;;; To run this program, download and install Chicken Scheme: https://call-cc.org
;;;;; Install necessary dependencies (as root):
;;;;;     chicken-install amb:2.1.6 loop:1.4 matchable:3.3 linear-algebra:1.4
;;;;; And then run it like this:
;;;;;     csi -ss amb.scm rho1 rho2 rho3
;;;;; For example (as in article):
;;;;;     csi -ss amb.scm 0.3333 0.3333 0.3333
;;;;; The program will then output the X matrix and the values of quality,
;;;;; finances and time.
(use srfi-1 amb amb-extras loop matchable linear-algebra)

;;;; misc;{{{
(define (partial f . args)
  (lambda rest
    (apply f (append args rest))))

(define (mapcat . args)
  (concatenate (apply map args)))

(define all-equal?
  (match-lambda
   [() #t]
   [(_) #t]
   [(x . ys) (every (lambda (y) (= x y)) ys)]
   [_ (error "An argument must be a list")]))

(define (sum lst)
  (apply + lst))

(define (sort-by key lst #!optional (less <))
  (sort lst (lambda (a b) (less (key a) (key b)))))
;}}}
;;;; graph;{{{
(define *graph-eq?* eqv?)

(define (graph-vertices graph)
  (apply lset-adjoin *graph-eq?* '() (concatenate graph)))

(define (graph-next graph v)
  (map second (filter (match-lambda ((v1 v2) (*graph-eq?* v1 v))) graph)))

(define (graph-last? graph v)
  (null? (graph-next graph v)))

(define (graph-seconds graph)
  (map second graph))

(define (graph-is-start? graph vertice)
  (not (member vertice (graph-seconds graph) *graph-eq?*)))

(define (graph-starts graph)
  (filter (partial graph-is-start? graph) (graph-vertices graph)))
;}}}
;;;; paths;{{{
(define (paths-from graph start)
  (if (graph-last? graph start)
      (list (list start))
      (let* ((after-start (graph-next graph start))
             (paths-from-after-start (mapcat (partial paths-from graph) after-start)))
        (map (partial cons start) paths-from-after-start))))

(define (paths graph)
  (let ((starts (graph-starts graph)))
    (mapcat (partial paths-from graph) starts)))
;}}}
;;;; matrix;{{{
(define (mat-sum mat)
  (let ((m (matrix-rows mat))
        (n (matrix-columns mat)))
    (loop for i from 0 to (- m 1) sum
          (loop for j from 0 to (- n 1) sum
                (matrix-ref mat i j)))))

(define (m*.+ m1 m2)
  (mat-sum (m*. m1 m2)))
;}}}
;;;; amb;{{{
(define (amb-bool)
  (amb 0 1))

(define (amb-bool-matrix m n)
  (list->matrix
   (list-tabulate m (lambda (i)
                      (list-tabulate n (lambda (j) (amb-bool)))))))

(define (amb-1-list n)
  (let ((lists (list-tabulate n (lambda (i) (append (make-list i 0)
                                               '(1)
                                               (make-list (- n i 1) 0))))))
    (choose lists)))
;}}}
;;;; optimization;{{{
(define (optimize f lst op)
  (car (sort-by f lst op)))

(define (best f lst)
  (optimize f lst >))

(define (worst f lst)
  (optimize f lst <))
;}}}
;;;; objective functions;{{{
(define (time* t P x)
  ;; changed time to time* to evade conflict with builtin
  (let ((contractors-count (matrix-columns t)))
    (max* (loop for p in P collect
                (loop for i in p sum
                      (loop for j from 0 to (- contractors-count 1) sum
                            (* (matrix-ref x (- i 1) j)
                               (matrix-ref t (- i 1) j))))))))

(define (quality q w x)
  (loop for i from 0 to (- (matrix-rows x) 1) sum
        (loop for j from 0 to (- (matrix-columns x) 1) sum
              (* (matrix-ref x i j)
                 (list-ref w i)
                 (matrix-ref q i j)))))

(define (finances f x)
  (m*.+ x f))
;}}}
;;;; input;{{{
(define oo 1e100) ; infinity

(define t
  (list->matrix `((5  3   ,oo)
                  (4  5   ,oo)
                  (4  3   ,oo)
                  (15 10  8  )
                  (15 6   ,oo)
                  (20 ,oo 3  )
                  (25 ,oo 4  )
                  (30 15  ,oo))))

(define f
  (list->matrix `((5  10  ,oo)
                  (5  5   ,oo)
                  (5  5   ,oo)
                  (10 5   5)
                  (25 15  ,oo)
                  (30 ,oo 10)
                  (40 ,oo 30)
                  (50 15  ,oo))))

(define q
  (list->matrix '((10 3 0)
                  (10 4 0)
                  (10 4 0)
                  (10 6 2)
                  (9  9 0)
                  (5  0 8)
                  (6  0 7)
                  (8  7 0))))

(define w '(30 20 10 20 10 5 3 2))

(define graph
  '((6 2)
    (5 2)
    (2 7)
    (2 8)
    (7 3)
    (8 3)
    (3 4)
    (4 1)))
;}}}
;;;; input-related functions;{{{
;;; Solve my problem (as defined in input section)
(define (make-w->min f min max)
  (lambda (x)
    (/ (- (f x) min)
       (- max min))))

(define (make-w->max f min max)
  (lambda (x)
    (/ (- max (f x))
       (- max min))))

(define (max* lst)
  (apply max lst))

(define (amb-make-x rows cols)
  (let ((lists (list-tabulate rows (lambda (_) (amb-1-list cols)))))
   (list->matrix lists)))

;;; generate all possible X matrices with given number of rows and columns
(define (make-all-xs rows cols)
  (amb-collect (amb-make-x rows cols)))

(define (my-solve rho1 rho2 rho3)
  (let* ((rows (matrix-rows t))
         (cols (matrix-columns t))
         (all-xs (make-all-xs rows cols))

         (quality (partial quality q w))
         (Qmin (quality (worst quality all-xs)))
         (Qmax (quality (best quality all-xs)))

         (time* (partial time* t (paths graph)))
         (Tmin (time* (worst time* all-xs)))
         (Tmax (time* (best time* all-xs)))

         (finances (partial finances f))
         (Fmin (finances (worst finances all-xs)))
         (Fmax (finances (best finances all-xs)))

         (w1 (make-w->max quality Qmin Qmax))
         (w2 (make-w->min finances Fmin Fmax))
         (w3 (make-w->min time* Tmin Tmax))

         (x (worst (lambda (x) (max (* rho1 (w1 x))
                                    (* rho2 (w2 x))
                                    (* rho3 (w3 x))))
                   all-xs))
         )
    `((x ,x)
      (q ,(quality x))
      (t ,(time* x))
      (f ,(finances x)))))
;}}}
;;; cli;{{{
(define (die-with-usage)
  (print "Usage: " (program-name) " rho1 rho2 rho3
where rho1, rho2 and rho3 are numbers")
  (exit 1))

(define (print-matrix m)
  (for-each print (map vector->list m)))

(define (main argv)
  (if (not (= (length argv) 3))
      (die-with-usage))
  (let ((rho1 (string->number (first  argv)))
        (rho2 (string->number (second argv)))
        (rho3 (string->number (third argv))))
    (if (and rho1 rho2 rho3)
        (let ((solution (my-solve rho1 rho2 rho3)))
          (print-matrix (alist-ref 'x solution))
          (newline)
          (print "Quality: " (alist-ref 'q solution))
          (newline)
          (print "Finances: " (alist-ref 'f solution))
          (newline)
          (print "Time: " (alist-ref 't solution))
          (newline))
        (die-with-usage))))
;}}}
; vim:foldmethod=marker
