(define (readthem)
  (do ((entry (read) (read))
       (l '() (cons entry l)))
      ((eof-object? entry) l)))

(define (uniquify l)
  (let loop ((l l) (done '()) (result '()))
    (if (null? l)
	result
	(let ((n (assoc 'name (car l))))
	  (if (or (not n) (member (cadr n) done))
	      (loop (cdr l) done result)
	      (loop (cdr l) (cons (cadr n) done) (cons (car l) result)))))))

(define (group-by-package l)
  (let loop ((l l) (groups '()))
    (if (null? l)
	groups
	(let* ((c (car l))
	       (p (assoc 'package c)))
	  (if p
	      (let ((e (assoc (cadr p) groups)))
		(if e
		    (begin (set-cdr! e (cons c (cdr e)))
			   (loop (cdr l) groups))
		    (loop (cdr l) (cons (cons (cadr p) (list c)) groups))))
	      ;; no group:
	      (loop (cdr l) (cons (cons #f (list c)) groups)))))))
	      
		      
	      


(define x (readthem))
(define y (group-by-package (uniquify x)))
