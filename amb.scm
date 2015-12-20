(use srfi-1 matchable loop)

;; misc
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

;; graph
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

;; paths
(define (paths-from graph start)
  (if (graph-last? graph start)
      (list (list start))
      (let* ((after-start (graph-next graph start))
             (paths-from-after-start (mapcat (partial paths-from graph) after-start)))
        (map (partial cons start) paths-from-after-start))))

(define (paths graph)
  (let ((starts (graph-starts graph)))
    (mapcat (partial paths-from graph) starts)))

;; matrix
(use linear-algebra)

(define (mat-sum mat)
  (let ((m (matrix-rows mat))
        (n (matrix-columns mat)))
    (loop for i from 0 to (- m 1) sum
          (loop for j from 0 to (- n 1) sum
                (matrix-ref mat i j)))))

(define (m*.+ m1 m2)
  (mat-sum (m*. m1 m2)))

;; amb
(use amb amb-extras)

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

;; asserts
(define (assert-list lst #!optional msg)
  (assert (list? lst) (or msg "Not a list")))

(define (assert-matrix m #!optional msg)
  (assert (matrix? m) (or msg "Not a matrix")))

(define (assert-same-dimensions m1 m2)
  (assert-matrix m1)
  (assert-matrix m2)
  (assert (and (= (matrix-columns m1) (matrix-columns m2))
               (= (matrix-rows m1) (matrix-rows m2)))))

;; ra
(define (ra #!key w t f q F T graph)
  (assert-list w)
  (assert-same-dimensions t f)
  (assert-same-dimensions t q)
  (assert (= (length w) (matrix-rows t)))
  (assert (> F 0))
  (assert (> T 0))
  (assert (list? graph))
  (let* ((bp-count (matrix-rows t))
         (contractors-count (matrix-columns t))
         (P (paths graph))
         (lists (list-tabulate bp-count (lambda (_) (amb-1-list contractors-count))))
         (x (list->matrix lists)))
    ;(assert-same-dimensions t x)
    ;; sum equals 1
    ;(loop for i from 0 to (- bp-count 1) do (required (= 1 (loop for j from 0 to (- contractors-count 1) sum (matrix-ref x i j)))))
    ;; time
    (loop for p in P do
          (required (<= (loop for i in p sum
                              (loop for j from 0 to (- contractors-count 1) sum
                                    (* (matrix-ref x (- i 1) j)
                                       (matrix-ref t (- i 1) j))))
                        T)))
    ;; finances
    (required (<= (m*.+ x f) F))
    x))

(define oo 1e100) ;; infinity

(define t
  (list->matrix `((5  3   ,oo ,oo ,oo)
                  (4  5   ,oo ,oo ,oo)
                  (4  3   ,oo ,oo ,oo)
                  (10 8   3   5   ,oo)
                  (15 6   ,oo 5   ,oo)
                  (20 ,oo 3   5   ,oo)
                  (20 ,oo 4   7   ,oo)
                  (30 15  ,oo 15  20))))

(define f
  (list->matrix `((5000  10000 ,oo   ,oo   ,oo)
                  (5000  5000  ,oo   ,oo   ,oo)
                  (5000  5000  ,oo   ,oo   ,oo)
                  (20000 15000 5000  10000 ,oo)
                  (25000 15000 ,oo   15000 ,oo)
                  (30000 ,oo   10000 15000 ,oo)
                  (40000 ,oo   30000 25000 ,oo)
                  (50000 15000 ,oo   15000 25000))))

(define q
  (list->matrix '((1   0.3 0   0   0)
                  (1   0.4 0   0   0)
                  (1   0.4 0   0   0)
                  (1   0.6 0.2 0.5 0)
                  (0.9 0.9 0   0.7 0)
                  (0.5 0   0.8 0.9 0)
                  (0.6 0   0.7 0.8 0)
                  (0.8 0.7 0   0.7 0.9))))

(define w '(0.3 0.2 0.1 0.2 0.1 0.05 0.03 0.02))

(define graph
  '((6 2)
    (5 2)
    (2 7)
    (2 8)
    (7 3)
    (8 3)
    (3 4)
    (4 1)))

(define (quality q w x)
  (loop for i from 0 to (- (matrix-rows x) 1) sum
        (loop for j from 0 to (- (matrix-columns x) 1) sum
              (* (matrix-ref x i j)
                 (list-ref w i)
                 (matrix-ref q i j)))))

(define (my-ra F T)
  (ra w: w t: t f: f q: q F: F T: T graph: graph))

(define (my-ra-best F T)
  (let ((ans (amb-collect (my-ra F T))))
    (car (sort-by (partial quality q w) ans >))))

;; cli
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
