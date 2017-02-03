


;; extracted from mat.lisp 1/17/2016 RJF

(defun backward ()
  (declare(special ax delta m rank))
  (do ((i (1- rank) (1- i)))
      ((< i 1))
    (do ((l (1+ rank) (1+ l)))
	((> l m))
      (setf (aref ax (aref *row* i) (aref *col* l))
	(let ((mess1  (pdifference
			 (ptimes (aref ax (aref *row* i) (aref *col* l))
				 (aref ax (aref *row* rank) (aref *col* rank)))
			 (do ((j (1+ i) (1+ j)) (sum 0))
			     ((> j rank) sum)
			   (setq sum (pplus sum (ptimes (aref ax (aref *row* i) (aref *col* j))
							(aref ax (aref *row* j) (aref *col* l))))))) )
	      (mess2 (aref ax (aref *row* i) (aref *col* i))  ))
	  (cond ((equal 0 mess1) 0)
		((equal 0 mess2) 0)
		(t   ;;   (pquotient mess1 mess2) ; fixed by line below. RJF 1/12/2017

		 (car (ratreduce mess1 mess2))
		 )
		))))
    (do ((l (1+ i) (1+ l)))
	((> l rank))
      (setf (aref ax (aref *row* i) (aref *col* l)) 0)))
  ;; PUT DELTA INTO THE DIAGONAL MATRIX
  (setq delta (aref ax (aref *row* rank) (aref *col* rank)))
  (do ((i 1 (1+ i)))
      ((> i rank))
    (setf (aref ax (aref *row* i) (aref *col* i)) delta))  )



(defun forward (*cpivot)
  (declare(special nvar rank nrow ax delta m))
  (setq delta 1)		  ;DELTA HOLDS THE CURRENT DETERMINANT
  (do ((k 1 (1+ k))
       (nvar nvar)   ;PROTECTS AGAINST TEMPORARAY RESETS DONE IN PIVOT
       (m m))
      ((or (> k nrow) (> k nvar)))
    (cond ((pivot ax k *cpivot) (return nil)))
    ;; PIVOT IS T IF THERE IS NO MORE NON-ZERO ROW LEFT. THEN GET OUT OF THE LOOP
    (do ((i (1+ k) (1+ i)))
	((> i nrow))
      (do ((j (1+ k) (1+ j)))
	  ((> j m))
	(setf (aref ax (aref *row* i) (aref *col* j))
	  
	  #+ignore
	   (pquotient 
		      delta)
	   
	   (let ((mess1 (pdifference (ptimes (aref ax (aref *row* k) (aref *col* k))
					       (aref ax (aref *row* i) (aref *col* j)))
				       (ptimes (aref ax (aref *row* i) (aref *col* k))
					       (aref ax (aref *row* k) (aref *col* j))))))
	     (cond ((equal mess1 0) 0)
		   ((equal delta 0) 0)
		   (t (car (ratreduce mess1 delta)))))
	       
	       
	       )))
    (do ((i (1+ k) (1+ i)))
	((> i nrow))
      (setf (aref ax (aref *row* i) (aref *col* k)) 0))
    (setq delta (aref ax (aref *row* k) (aref *col* k))))
  ;; UNDOES COLUMN HACK IN PIVOT.
  (or *cpivot (do ((i 1 (1+ i))) ((> i m)) (setf (aref *col* i) i)))
  (setq rank (min nrow nvar)))

;; other places that pquotient is called 
#|
rootfac in sinint.lisp
resprog
pdecomp*  in combin.lisp

|#