;; sanitize .info file
(define verboten '("djbfft-0.76" "athfft" "pfftw"))

(do ((entry (read) (read)))
    ((eof-object? entry))
  (if (and (assoc 'name entry) (member (cadr (assoc 'name entry)) verboten))
      'nothing-to-do
      (begin   (display "(")
               (newline)
               (for-each (lambda (x) (write x) (newline)) entry)
               (display ")")
               (newline))))
