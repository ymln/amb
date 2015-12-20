(use srfi-1 matchable loop linear-algebra amb amb-extras)

;; misc;{{{
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
;; graph;{{{
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
;; paths;{{{
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
;; matrix;{{{
(define (mat-sum mat)
  (let ((m (matrix-rows mat))
        (n (matrix-columns mat)))
    (loop for i from 0 to (- m 1) sum
          (loop for j from 0 to (- n 1) sum
                (matrix-ref mat i j)))))

(define (m*.+ m1 m2)
  (mat-sum (m*. m1 m2)))
;}}}
;; amb;{{{
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
;; asserts;{{{
(define (assert-list lst #!optional msg)
  (assert (list? lst) (or msg "Not a list")))

(define (assert-matrix m #!optional msg)
  (assert (matrix? m) (or msg "Not a matrix")))

(define (assert-same-dimensions m1 m2)
  (assert-matrix m1)
  (assert-matrix m2)
  (assert (and (= (matrix-columns m1) (matrix-columns m2))
               (= (matrix-rows m1) (matrix-rows m2)))))

(define (max* lst)
  (apply max lst))

(define (amb-make-x rows cols)
  (let ((lists (list-tabulate rows (lambda (_) (amb-1-list cols)))))
   (list->matrix lists)))

(define (assert-sum-equals-1 matrix)
  (let ((rows (matrix-rows matrix))
        (cols (matrix-columns matrix)))
    (loop for i from 0 to (- rows 1) do
          (assert (= 1
                     (loop for j from 0 to (- cols 1) sum (matrix-ref matrix i j)))))))
(define (ra-asserts w t f q F T graph x)
  (assert-list w)
  (assert-same-dimensions t f)
  (assert-same-dimensions t q)
  (assert (= (length w) (matrix-rows t)))
  (assert (> F 0))
  (assert (> T 0))
  (assert (list? graph))
  (assert-sum-equals-1 x)
  (assert-same-dimensions t x))
;}}}
;; optimization;{{{
(define (optimize f lst op)
  (car (sort-by f lst op)))

(define (best f lst)
  (optimize f lst >))

(define (worst f lst)
  (optimize f lst <))
;}}}
;; objective functions;{{{
(define (time* x t P)
  ; changed time to time- to evade conflict with builtin
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

(define (finances x f)
  (m*.+ x f))
;}}}
;; ra;{{{
(define (ra #!key w t f q F T graph)
  (let* ((bp-count (matrix-rows t))
         (contractors-count (matrix-columns t))
         (x (amb-make-x bp-count contractors-count)))
    (ra-asserts w t f q F T graph x)
    (required (<= (time* x t (paths graph)) T))
    (required (<= (finances x f) F))
    x))
;}}}
;; input;{{{
(define oo 1e100) ;; infinity

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
;; input-related functions;{{{
(define (my-ra F T)
  (ra w: w t: t f: f q: q F: F T: T graph: graph))

(define (my-ra-best F T)
  (let ((ans (amb-collect (my-ra F T))))
    (best (partial quality q w) ans)))
;}}}
;; cli;{{{
(define (die-with-usage)
  (print "Usage: " (first (argv)) " F T
where F and T are numbers")
  (exit 1))

(define (main)
  (if (not (= (length (argv)) 4))
      (die-with-usage))
  (let ((F (string->number (third  (argv))))
        (T (string->number (fourth (argv)))))
    (if (and F T)
        (let ((best (my-ra-best F T)))
          (for-each print (vector->list best))
          (newline)
          (print "Quality: " (quality q w best)))
        (die-with-usage))))

(main)
(exit)
;}}}
; vim:foldmethod=marker
