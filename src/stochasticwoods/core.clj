(ns stochasticwoods.core
  (:use [clojure.repl])
  (:require [clojure.core.rrb-vector :as rrb]))

;; Decision Trees
(defn- square
  [n]
  (* n n))

(defn gini-impurity
  [category_seq]
  (let [n (count category_seq)
        f (frequencies category_seq)]
    (- 1 (reduce + (map (fn [[k v]] (square (/ v n))) f)))))

(defn categorical?
  [data]
  (if (every? number? data)
    false
    true))

(defn infer-datatypes
  [data]
  (let [ks (if (map? (first data)) (keys (first data))
               (range (count (first data))))]
    (apply merge
           (map #(if (categorical? (map (fn [row] (get row %)) data))
                   {% :categorical}
                   {% :numeric})
                ks))))

(defn infer-possible-vals
  ([data] (infer-possible-vals data
                               (if (map? (first data))
                                 (keys (first data))
                                 (range (count (first data))))))
  ([data ks]
   (let [f (fn [k] (apply hash-set (map #(get % k) data)))]
     (apply merge (map (fn [k] {k (f k)}) ks)))))

(defn- break-fns
  ([data ks] (break-fns data ks (infer-datatypes data)))
  ([data ks datatypes]
   (let [datasets (infer-possible-vals data ks)]
     (fn [n]
       (map (fn [p]
              (if (= :categorical (datatypes n))
                [(fn [x] (= (get x n) p)) (str "[" n "] == " p)]
                [(fn [x] (<= (get x n) p)) (str "[" n "] <= " p)]))
            (datasets n))))))

;; (defn- partition-data-by-breakfn
;;   "This ugly function sorts the data into
;;   two bins using the break function."
;;   [data response]
;;   (let [m (count data)]
;;     (fn [[bfn s]]
;;       (loop [leftd (list)
;;              rightd (list)
;;              leftr (list)
;;              rightr (list)
;;              i 0]
;;         (if (< i m)
;;           (if (bfn (nth data i))
;;             (recur (conj leftd (nth data i))
;;                    rightd
;;                    (conj leftr (nth response i))
;;                    rightr
;;                    (inc i))
;;             (recur leftd
;;                    (conj rightd (nth data i))
;;                    leftr
;;                    (conj rightr (nth response i))
;;                    (inc i)))
;;           [[leftd leftr] [rightd rightr]])))))

(defn- special-assoc!
  [m k v]
  (assoc! m k (inc (get m k 0))))

(defn- partition-data-by-breakfn
  "This ugly function sorts the data into
  two bins using the break function."
  [data response]
  (let [m (count data)]
    (fn [[bfn s]]
      (loop [leftd (list)
             rightd (list)
             leftr (list)
             rightr (list)
             lvals (transient {})
             rvals (transient {})
             i 0]
        (if (< i m)
          (let [d (nth data i)
                r (nth response i)]
            (if (bfn d)
              (recur (conj leftd d)
                     rightd
                     (conj leftr r)
                     rightr
                     (special-assoc! lvals r 1)
                     rvals
                     (inc i))
              (recur leftd
                     (conj rightd d)
                     leftr
                     (conj rightr r)
                     lvals
                     (special-assoc! rvals r 1)
                     (inc i))))
          [[leftd leftr (persistent! lvals)]
           [rightd rightr (persistent! rvals)]])))))

;; (defn gini-impurity-after-break
;;   [[[ld lr] [rd rr]]]
;;   ;;(* 0.5 (+  (gini-impurity lr) (gini-impurity rr)))
;;   ;;(max (gini-impurity lr) (gini-impurity rr))
;;   (let [n (+ (count lr) (count rr))]
;;     (+ (* (gini-impurity lr) (/ (count lr) n))
;;        (* (gini-impurity rr) (/ (count rr) n)))))

(defn- gifromfreqs
  [f]
  (let [n (reduce + (map second f))]
    (- 1 (reduce + (map (fn [[k v]] (square (/ v n))) f)))))

(defn gini-impurity-after-break
  [[[ld lr lv] [rd rr rv]]]
  ;;(* 0.5 (+  (gini-impurity lr) (gini-impurity rr)))
  ;;(max (gini-impurity lr) (gini-impurity rr))
  (let [n (+ (count lr) (count rr))]
    (+ (* (gifromfreqs lv)
          (/ (count lr) n))
       (* (gifromfreqs rv)
          (/ (count rr) n)))))

(defn grab-details
  [data response breakfns]
  (map #(let [partition_fn (partition-data-by-breakfn data response)
              p (partition_fn %)]
          (list % p (gini-impurity-after-break p)))
       breakfns))

(defn build-decision-tree
  ([data response]
   (build-decision-tree data
                        response
                        (if (map? (first data))
                          (keys (first data))
                          (range (count (first data))))
                        (infer-datatypes data)))
  ([data response columnkeys datatypes]
   (let [gi0 (gini-impurity response)
         m (count data)
         n (count columnkeys)
         foo (break-fns data columnkeys datatypes)
         breakfns (mapcat foo columnkeys)
         allbreaks (grab-details data response breakfns)
         [mf [[lpd lpr] [rpd rpr]] gi] (reduce (partial min-key last)
                                               allbreaks)]
     (cond
       (<= gi0 gi) (frequencies response)
       :else (list mf
                   (build-decision-tree (into [] lpd) (into [] lpr)
                                        columnkeys datatypes)
                   (build-decision-tree (into [] rpd) (into [] rpr)
                                        columnkeys datatypes))))))

;; Printing trees
(defn adjacency-list
  ([tree] (adjacency-list tree 1))
  ([tree n]
   (if (not (list? tree)) [n tree]
       (loop [[s1 & sr] (rest tree)
              ret []
              m n]
         (if (nil? s1) ret
             (let [v (adjacency-list s1 (inc m))]
               (if (number? (first v))
                 (recur sr
                        (conj ret [[n (first v)] (second v)])
                        (first v))
                 (recur sr
                        (concat (conj ret
                                      [[n (first (first (first v)))]
                                       (first s1)]) v)
                        (apply max (mapcat first v))))))))))

(defn dot-format
  [tree]
  (let [al (adjacency-list tree)
        edges (map first al)
        fs (fn [id node]
             (str id
                  " [label=\""
                  (if (map? node) node (second node))
                  "\", shape=\"box\"];\n"))]
    (str
     "digraph Tree {\n"
     (fs 1 (first tree))
     (reduce str
             (map (fn [[[e1 e2] node]] (fs e2 node)) al))
     (reduce str
             (map (fn [[e1 e2]] (str e1 " -> " e2 ";\n"))
                  edges))
     "}")))

;; Predictions
(defn predict-from-tree
  [[[rootfn rootlabel] & subtree] x]
  (let [nxt {true 0 false 1}
        newtree (nth subtree (nxt (rootfn x)))]
    (if (map? newtree) newtree (recur newtree x))))

(defn predict-from-forest
  [forest x]
  (reduce (partial merge-with +)
          (map #(predict-from-tree % x) forest)))

(defn classify-from-tree
  [tree x]
  (first (reduce (partial max-key second)
                 (predict-from-tree tree x))))

(defn classify-from-forest
  [forest x]
  (first (reduce (partial max-key second)
                 (predict-from-forest forest x))))

;; Evaluation
(defn train-and-validate
  [traindata trainresponse testdata testresponse]
  (let [dtree (build-decision-tree traindata trainresponse)
        [tpos fpos tneg fneg] [0 0 0 0]
        f (fn [x r] (if (= (classify-from-tree dtree x) r)
                     :correct
                     :incorrect))]
    (frequencies (map f testdata testresponse))))

(defn- drop-from-rrb-vector
  [v index]
  (rrb/catvec (rrb/subvec v 0 index)
              (rrb/subvec v (inc index))))

(defn sample-with-replacement
  [s n]
  (let [v (if (vector? s) s (into [] s))
        m (count v)]
    (repeatedly n #(nth v (rand-int m)))))

(defn sample-without-replacement
  [things num]
  (let [v (if (vector? things) things (into [] things))
        m (min num (count v))]
    (loop [curr v
           n m
           ret (list)]
      (cond
        (zero? n) ret
        :else (let [i (rand-int (count curr))]
                (recur (drop-from-rrb-vector curr i)
                       (dec n)
                       (conj ret (curr i))))))))

(defn build-forest
  [data response]
  (let [ks (vec (if (map? (first data)) (keys (first data))
                    (range (count (first data)))))
        ntrees 500
        datatypes (infer-datatypes data)
        indices (vec (range (count data)))
        nsamples (int (/ (count data) 1))
        nfeatures (min (count ks)
                       (+ 2 (int (java.lang.Math/sqrt (count ks)))))
        samplefn (fn [] [(sample-with-replacement indices nsamples)
                        (sample-without-replacement ks nfeatures)])
        datafn (fn [sinds] (mapv (partial get data) sinds))
        responsefn (fn [sinds] (mapv (partial get response) sinds))
        f (fn [[sind features]]
            (build-decision-tree (datafn sind)
                                 (responsefn sind)
                                 features
                                 datatypes))]
    (take ntrees
          (pmap #(hash-map :features (second %) :tree (f %))
                (repeatedly samplefn)))))

(defn train-and-validate
  [traindata trainresponse testdata testresponse]
  (let [forest (build-forest traindata trainresponse)
        [tpos fpos tneg fneg] [0 0 0 0]
        f (fn [x r] (if (= (classify-from-forest (map :tree forest) x) r)
                     :correct
                     :incorrect))]
    (frequencies (map f testdata testresponse))))

(time (train-and-validate (into [] (take 300 data))
                          (into [] (take 300 response))
                          (into [] (nthrest data 300))
                          (into [] (nthrest response 300))))


(def forest (build-forest data response))
(time (count forest))


(def t (build-decision-tree data response))
(classify-from-tree t [:female 100])
(time (predict-from-forest (repeat 1000 t) [:male 50]))

(def blah (read-string (slurp "/home/kanderson/Downloads/test.edn")))
(def data (map #(dissoc % :new_tumor_event_after_initial_treatment)
               blah))
(def response (map :new_tumor_event_after_initial_treatment blah))

(def f "/home/kanderson/Code/Bioinformatics/lung_risk/missing_filled_with_year.csv")
(def foo (with-open [rdr (clojure.java.io/reader f)]
           (let [l (line-seq rdr)
                 ks (map keyword
                         (clojure.string/split (first l)
                                               #","))
                 f (fn [row] (map read-string
                                 (clojure.string/split row #",")))
                 f2 (fn [row] (zipmap ks (f row)))]
             (mapv f2 (rest l)))))
(def data (mapv #(dissoc % :new_tumor_event_after_initial_treatment)
                foo))
(def response (mapv :new_tumor_event_after_initial_treatment foo))
(time (build-decision-tree data response))

(def data
  (into []
        (repeatedly 3000 #(if (> 0.5 (rand))
                            [:male (rand-int 100) (rand) (rand)]
                            [:female (rand-int 100) (rand) (rand)]))))
(def response
  (mapv (fn [[s a]] (if (= s :male)
                     (if (> a 80)
                       :yes
                       (if (< (rand) 0.8)
                         :no
                         :yes))
                     :no))
        data))

(infer-datatypes data)
(frequencies response)
(time (build-decision-tree data response))
(second (adjacency-list (build-decision-tree data response)))
(spit "/home/kanderson/Downloads/foobar.dot"
      (dot-format (build-decision-tree data response
                                       :columnkeys [:tumor_weight])))



