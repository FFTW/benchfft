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
	      

(define (every? pred? l)
  (or (null? l) 
      (and (pred? (car l)) (every? pred? (cdr l)))))

(define (common? info-entry package-members)
  (every? (lambda (m) (member info-entry m)) package-members))

(define (filter pred? l)
  (if (null? l)
      l
      (if (pred? (car l))
	  (cons (car l) (filter pred? (cdr l)))
	  (filter pred? (cdr l)))))


(define (do-package package)
  (let* ((name (car package))
	 (package-members (cdr package))
	 (common-entries
	  (filter (lambda (info-entry) 
		    (common? info-entry package-members))
		  (car package-members)))
	 (specific-entries
	  (map (lambda (m) 
		 (filter (lambda (info-entry) 
			   (not (common? info-entry package-members)))
			 m))
	       package-members)))
    (htmlize name common-entries specific-entries)))

(define packages (group-by-package (uniquify (readthem))))

(define (htmlize name common-entries specific-entries)
  'todo)

(for-each do-package packages)