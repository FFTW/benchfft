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

(define (name2package n)
  (if n
      (let* ((ns (cadr n))
	     (hyphen (string-index ns #\-)))
	(list 'package (if hyphen
			   (substring ns 0 hyphen)
			   ns)))
      #f))

(define (group-by-package l)
  (let loop ((l l) (groups '()))
    (if (null? l)
	groups
	(let* ((c (car l))
	       (p0 (assoc 'package c))
	       (p (if p0 p0 (name2package (assoc 'name c)))))
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

(define (writeln . l)
  (for-each display l)
  (newline))

(define (htmlize package-name common-entries specific-entries)
  (begin (htmlize-one package-name common-entries)
	 (if (not (every? null? specific-entries))
	     (begin
	       (writeln "<ul>")
	       (for-each (lambda (e) (if (not (null? e)) (htmlize-one #f e)))
			 specific-entries)
	       (writeln "</ul>")))))

(define (assoc* sym l)
  (filter (lambda (x) (eq? sym (car x))) l))

(define (htmlize-one name entries)
  (let ((name (or name 
		  (and (assoc 'name entries) (cadr (assoc 'name entries)))
		  "unknown")))
    (writeln "<li> " name)
    (writeln "<ul>")
    (let ((url (assoc 'url entries))
	  (url-was-valid-on (assoc 'url-was-valid-on entries))
	  )
      (if url
	  (begin
	    (writeln "<li>URL: <a href=\"" (cadr url) "\">" (cadr url) "</a>")
	    (if url-was-valid-on
		(writeln "(was valid on " (cadr url-was-valid-on) ")"))))
      (maybe-plural "Author" "Authors" 'author entries)
      (maybe "Year" 'year entries)
      (maybe "Version" 'version entries)
      (maybe-plural "Language" "Languages" 'language entries)
      (maybe "References" 'bibitem entries)
      (for-each (lambda (note) (writeln "<li>Note: " (cadr note)))
		(assoc* 'notes entries))
      )
    (writeln "</ul>") ))

(define (maybe name sym entries)
  (let ((x (assoc sym entries)))
    (if x (writeln "<li>" name ": " (cadr x)))))

(define (maybe-plural singular plural sym entries)
  (let ((things (assoc* sym entries)))
    (cond ((= (length things) 1) 
	   (apply writeln "<li>" singular ": " (map cadr things)))
	  ((> (length things) 1)
	   (apply writeln "<li>" plural ": " 
		  (cadr (car things))
		  (map (lambda (x) (string-append ", " (cadr x)))
		       (cdr things)))))))

(define (compare-packages p1 p2)
  (string-ci<? (string-downcase (car p1)) (string-downcase (car p2))))

(writeln "<html><body>")
(writeln "<ul>")
(for-each do-package (sort packages compare-packages))
(writeln "</ul>")
(writeln "</body></html>")


