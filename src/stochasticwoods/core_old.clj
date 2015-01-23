(ns stochasticwoods.core
  (:use [clojure.repl]))

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

;; (defn- break-fns
;;   ([data ks] (break-fns data ks identity (infer-datatypes data)))
;;   ([data ks labels] (break-fns data ks labels (infer-datatypes data)))
;;   ([data ks labels datatypes]
;;    (let [datasets (infer-possible-vals data ks)]
;;      (fn [n]
;;        (map (fn [p]
;;               (if (= :categorical (datatypes n))
;;                 [(fn [x] (= (get x n) p)) (str "[" (labels n) "] == " p) (labels n)]
;;                 [(fn [x] (<= (get x n) p)) (str "[" (labels n) "] <= " p) (labels n)]))
;;             (datasets n))))))


(defn- break-fns
  ([data ks] (break-fns data ks (infer-datatypes data)))
  ([data ks datatypes]
   (let [datasets (infer-possible-vals data ks)]
     (fn [n]
       (map (fn [p]
              (if (= :categorical (datatypes n))
                (list (fn [x] (= (get x n) p)) n p)
                (list (fn [x] (<= (get x n) p)) n p)))
            (datasets n))))))

(defn- partition-data-by-breakfn
  "This ugly function sorts the data into
  two bins using the break function."
  [data response]
  (let [m (count data)]
    (fn [[bfn]]
      (loop [leftd (list)
             rightd (list)
             leftr (list)
             rightr (list)
             lvals (transient {})
             rvals (transient {})
             i (int 0)]
        (if (< i m)
          (let [d (nth data i)
                r (nth response i)]
            (if (bfn d)
              (recur (conj leftd d)
                     rightd
                     (conj leftr r)
                     rightr
                     (assoc! lvals r (inc (get lvals r (int 0))))
                     rvals
                     (inc i))
              (recur leftd
                     (conj rightd d)
                     leftr
                     (conj rightr r)
                     lvals
                     (assoc! rvals r (inc (get rvals r (int 0))))
                     (inc i))))
          [[leftd leftr (persistent! lvals)]
           [rightd rightr (persistent! rvals)]])))))

(defn- gifromfreqs
  [f ^Double n]
  (- 1
     (apply +
            (map (fn [[k ^Double v]]
                   (let [vv (/ v n)]
                     (* vv vv)))
                 f))))

(defn- gini-impurity-after-break
  [[[ld lr lv] [rd rr rv]]]
  (let [nlr (double (count lr))
        nrr (double (count rr))
        n (double (+ nlr nrr))]
    (+ (* (gifromfreqs lv nlr)
          (/ nlr n))
       (* (gifromfreqs rv nrr)
          (/ nrr n)))))

(defn- grab-details
  [data response breakfns]
  (map #(let [partition_fn (partition-data-by-breakfn data response)
              p (partition_fn %)]
          (list % p (gini-impurity-after-break p)))
       breakfns))

(defn- build-split-fn
  [datatypes labels]
  (fn [[f cat_index split]]
    (let [l (labels cat_index)]
      {:function (if (= :categorical (datatypes cat_index))
                   `(fn [~'x] (= (get ~'x ~cat_index) ~split))
                   `(fn [~'x] (<= (get ~'x ~cat_index) ~split)))
       :label (if (= :categorical (datatypes cat_index))
                (str l " == " split)
                (str l " <= " split))
       :variable l})))

(defn build-decision-tree
  ([data response]
   (build-decision-tree data response
                        (if (map? (first data))
                          (keys (first data))
                          (range (count (first data))))))
  ([data response columnkeys]
   (build-decision-tree data response columnkeys identity))
  ([data response columnkeys labels]
   (build-decision-tree data response columnkeys labels
                        (infer-datatypes data)))
  ([data response columnkeys labels datatypes]
   (let [gi0 (double (gini-impurity response))
         m (count data)
         n (count columnkeys)
         bfn (break-fns data columnkeys datatypes)
         breakfns (mapcat bfn columnkeys)
         allbreaks (grab-details data response breakfns)
         splitfn (build-split-fn datatypes labels)
         [mf
          [[lpd lpr lv]
           [rpd rpr rv]]
          gi :as d] (reduce (partial min-key last)
                            allbreaks)]
     (cond
       (< (java.lang.Math/abs
           (- gi0 gi)) 1E-5) {:freqs (frequencies response)
                              :gini-impurity gi0}
           :else (list {:split-fn (splitfn mf)
                        :gini-impurity gi0
                        :freqs (merge-with + lv rv)}
                       (build-decision-tree (vec lpd) (vec lpr)
                                            columnkeys labels datatypes)
                       (build-decision-tree (vec rpd) (vec rpr)
                                            columnkeys labels datatypes))))))

;; Printing trees
(defn adjacency-list
  ([tree] (if (list? tree)
            (adjacency-list tree 1)))
  ([tree n]
   (if (not (list? tree)) [n tree]
       (let [[root left right] tree
             lal (adjacency-list left (inc n))
             maxfn (fn [al] (if (number? (first al)) (first al)
                               (apply max (mapcat first al))))
             m (maxfn lal)
             ral (adjacency-list right (inc m))
             f (fn [t al side] (if (number? (first al))
                                [[[n (first al)] (second al) side]]
                                (conj al
                                      [[n (first (first (first al)))]
                                       (first t) side])))]
         (concat (f left lal true) (f right ral false))))))

(defn dot-format
  ([tree]
   (let [al (adjacency-list tree)
         edges (map first al)
         sides (map #(nth % 2) al)
         fs (fn [id node]
              (str id
                   " [label=\""
                   (if (contains? node :split-fn)
                     (str (get-in node [:split-fn :label]) ", "))
                   (node :freqs)
                   (if (contains? node :split-fn)
                     "\",shape=\"box")
                   "\"];\n"))]
     (str
      "digraph Tree {\n"
      (fs 1 (first tree))
      (reduce str
              (map (fn [[[e1 e2] node]] (fs e2 node)) al))
      (reduce str
              (map (fn [[e1 e2] s] (str e1 " -> " e2
                                       " [label=\""
                                       s
                                       "\"];\n"))
                   edges sides))
      "}"))))

(defn show-tree
  [tree]
  (let [dot (dot-format tree)
        pdf (:out (clojure.java.shell/sh
                   "dot" "-Tpdf" :in dot :out-enc :bytes))]
    (clojure.java.shell/sh "display" :in pdf)))

;; Predictions
(defn predict-from-tree
  [tree x]
  (if (map? tree) (:freqs tree)
      (let [[root left right] tree
            nxt {true 0 false 1}
            f (get root :internal_fn)
            result (if (nil? f)
                     ((eval (get-in root[:split-fn :function])) x)
                     (f x))]
        (if result
          (recur left x)
          (recur right x)))))

(defn predict-from-forest
  [forest x]
  (reduce (partial merge-with +)
          (map #(predict-from-tree (:tree %) x) forest)))

(defn classify-from-tree
  [tree x]
  (first (reduce (partial max-key second)
                 (predict-from-tree tree x))))

(defn classify-from-forest
  [forest x]
  (first (reduce (partial max-key second)
                 (predict-from-forest forest x))))

;; Evaluation
(defn sample-with-replacement
  [s n]
  (let [v (if (vector? s) s (into [] s))
        m (count v)]
    (repeatedly n #(nth v (rand-int m)))))

(defn sample-without-replacement
  [things n]
  (take n (shuffle things)))

(defn build-forest
  [data response labels]
  (let [ks (vec (if (map? (first data)) (keys (first data))
                    (range (count (first data)))))
        ntrees 500
        datatypes (infer-datatypes data)
        indices (vec (range (count data)))
        nsamples (int (/ (count data) 1))
        nfeatures (min (count ks)
                       (+ 4 (int (java.lang.Math/sqrt (count ks)))))
        samplefn (fn [] [(sample-with-replacement indices nsamples)
                        (sample-without-replacement ks nfeatures)])
        dfn (fn [d sinds] (mapv (partial get d) sinds))
        f (fn [[sind features]]
            (build-decision-tree (dfn data sind)
                                 (dfn response sind)
                                 features
                                 (zipmap features (dfn labels features))
                                 datatypes))]
    (take ntrees
          (map #(do (println %2)
                     (hash-map :features (second %1)
                              :indices (first %1)
                              :tree (f %1)))
               (repeatedly samplefn)
               (range)))))


(defn roc-curve
  [forest category_key testdata testresponses]
  (let [f (fn [m] (/ (get m category_key 0)
                    (apply + (map second m))))
        l (map (comp f (partial predict-from-forest forest))
               testdata)
        f3 (fn [pred actual] (= pred actual))
        f4 (fn [pred actual] (if (= pred category_key)
                              (if (= pred actual)
                                {:true-positive 1}
                                {:false-positive 1})
                              (if (= actual category_key)
                                {:false-negative 1}
                                {:true-negative 1})))
        f2 (fn [threshold] (apply merge-with +
                                 (map #(f4 (if (>= %1 threshold)
                                             category_key
                                             nil)
                                           %2)
                                      l testresponses)))
        errors (map f2 (range 0 1 0.1))
        n (count testresponses)
        tot+ ((frequencies
               (map #(= % category_key) testresponses))
              true)
        tot- (- n tot+)]
    [(map #(/ (get % :false-positive 0) tot-) errors)
     (map #(/ (get % :true-positive 0) tot+) errors)]))

(defn calc-auc
  [[false_pos true_pos]]
  (apply +
         (map (fn [[fp1 fp2] [tp1 tp2]]
                (* (- fp1 fp2)
                   (* 0.5 (+ tp1 tp2))))
              (partition 2 1 false_pos)
              (partition 2 1 true_pos))))

(defn train-and-validate
  [data response train_indices test_indices labels category_key]
  (let [grab (fn [indices data] (mapv #(get data %) indices))
        traind (grab train_indices data)
        trainr (grab train_indices response)
        testd (grab test_indices data)
        testr (grab test_indices response)
        forest (build-forest traind trainr labels)]
    (roc-curve forest category_key testd testr)))

(defn cross-validate
  [data response fold labels category_key]
  (let [n (count response)
        m (int (/ n fold))
        f (map (fn [x] [(concat (range 0 (* m x))
                               (range (* m (inc x)) n))
                       (range (* m x) (* m (inc x)))])
               (range fold))]
    (map (fn [[train test]]
           (train-and-validate data response train
                               test labels category_key))
         f)))

(defn permute-variable
  [data feature_key]
  (let [v (map #(get % feature_key) data)
        f (fn [row pval] (assoc row feature_key pval))]
    (mapv f data (shuffle v))))

(defn- var-imp
  [data tree oob_indices]
  (let [f (fn [d] (map (fn [i] (classify-from-tree tree (get d i)))
                      oob_indices))
        c (f data)
        f2 (fn [c1 c2] (= c1 c2))]
    (fn [feature_key]
      (let [m
            (apply merge-with
                   + (repeatedly 20
                                 #(frequencies
                                   (map f2 c
                                        (f (permute-variable
                                            data
                                            feature_key))))))]
        {feature_key (float
                      (/ (get m false 0)
                         (+ (get m false 0) (get m true 0))))}))))

(defn- reify-fns-in-tree
  [tree]
  (if (map? tree) tree
      (let [[root left right] tree
            fun (get-in root [:split-fn :function])
            f (eval fun)]
        (list (assoc root :internal_fn f)
              (reify-fns-in-tree left)
              (reify-fns-in-tree right)))))

(defn- var-imp-tree
  [data]
  (let [all (apply hash-set (range (count data)))]
    (fn
      [{tree :tree indices :indices features :features}]
      (let [t (reify-fns-in-tree tree)
            ib (apply hash-set indices)
            oob (clojure.set/difference all ib)
            f (var-imp data t oob)]
        (apply merge
               (map f features))))))

(defn gini-decrease
  [{t :tree}]
  (letfn [(func [tree]
            (if (not (list? tree))
              [(tree :gini-impurity) {}]
              (let [[root left right] tree
                    gi (root :gini-impurity)
                    lsubtree (func left)
                    rsubtree (func right)
                    mean (fn [a b] (* 0.5 (+ a b)))
                    f (fn [[g rem]]
                        (merge-with +
                                    {(get-in root [:split-fn :variable])
                                     (float (- gi g))}
                                    rem))]
                [gi (merge-with + (f lsubtree) (f rsubtree))])))]
    (second (func t))))


(defn calc-var-importance
  [data forest & {:keys [importance_fn]
                  :or {importance_fn (var-imp-tree data)}}]
  (let [mean (fn [a b] (* 0.5 (+ a b)))]
    (apply merge-with mean (map #(do (println %2) (importance_fn %1))
                                 forest (range)))))



(def var_importance (sort-by second (time (calc-var-importance data2 forest2))))
(def gini_importance (sort-by second (time (calc-var-importance data2 forest2 :importance_fn gini-decrease))))

(take 100 var_importance)
(take 100 gini_importance)

(float
 (- 1
    (/ (* 6
          (apply +
                 (map #(square (- (first %) (second %)))
                      (for [[i r1] (map list var_importance (range))]
                        (vector r1
                                (last (first
                                       (remove nil?
                                               (filter (fn [[j r2]] (= (first i) (first j)))
                                                       (map list gini_importance (range)))))))))))
       (* (count var_importance)
          (dec (square (count var_importance)))))))

(def gini_importance (into [] (concat [(first var_importance)] gini_importance)))
(first var_importance)
(println gini_importance)
(second
 (for [[i r1] (map list var_importance (range))]
   (filter (fn [[j r2]] (= (first i) (first j)))
           (map list gini_importance (range)))))

(doseq [a (sort-by second (time (calc-var-importance data2 forest2)))]
  (println a))

(def forest2 (build-forest data2 response2
                           (zipmap (keys (first data2))
                                   (keys (first data2)))))

(let [foo (into [] (sort-by second
                            (time (calc-var-importance data forest
                                                       :importance_fn gini-decrease))))
      n (count foo)]
  (doseq [a (subvec foo (- n 100))]
    ;;(println (gene_list (first a)) (second a))
    (println a)
    ))

[data response train_indices test_indices labels category_key]



(def forest (build-forest data
                          response
                          ;;(zipmap (range) gene_list)
                          (zipmap (keys (first data))
                                  (keys (first data)))))
(time (count forest))
(println (keys (first data)))

(defn special-read-string
  [x]
  (let [v (read-string x)]
    (if (symbol? v) (keyword v)
        v)))

(def f "/home/kanderson/Code/Bioinformatics/lung_risk/missing_filled_with_year.csv")
(def foo (with-open [rdr (clojure.java.io/reader f)]
           (let [l (line-seq rdr)
                 ks (map keyword
                         (clojure.string/split (first l)
                                               #","))
                 f (fn [row] (map special-read-string
                                 (clojure.string/split row #",")))
                 f2 (fn [row] (zipmap ks (f row)))]
             (mapv f2 (rest l)))))

(def data2 (mapv #(dissoc % :new_tumor_event_after_initial_treatment)
                foo))
(def response2 (mapv :new_tumor_event_after_initial_treatment foo))
(build-decision-tree data2 response2)
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


(def f "/home/kanderson/Code/Bioinformatics/lung_risk/bar.csv")
(def foo (with-open [rdr (clojure.java.io/reader f)]
           (let [l (line-seq rdr)
                 meta (clojure.string/split (first l) #", ")
                 f (fn [row] (mapv special-read-string
                                  (clojure.string/split row #", ")))]
             [meta (mapv f (rest l))])))
(defn mean [row]
  (let [n (count row)]
    (/ (reduce + row) n)))
(map mean
     data)
(time (dotimes [n 100000] ((fn [x y]
                             (let [a (long x)
                                   b (long y)]
                               (* a b))) 1 2)))
(def gene_list (into [] (rest (first foo))))
(def data (mapv #(into [] (rest %)) (second foo)))
(def response (mapv first (second foo)))

(def forest (build-forest data
                          response
                          (zipmap (range) gene_list)))
(time (count forest))
(show-tree ((nth forest2 4) :tree))


(use 'incanter.core)
(use 'incanter.charts)
(calc-auc (roc-curve forest2 (symbol "YES") data2 response2))
(let [[fp tp] (roc-curve forest2 (symbol "YES") data2 response2)]
  (view (xy-plot fp tp)))


(defn incanter-plot-roc
  [chart [false_positives true_positives]]
  (add-lines chart false_positives true_positives))

(view (reduce incanter-plot-roc
              (xy-plot (first (first foo))
                       (second (first foo)))
              (rest foo)))


(def foo (cross-validate data2 response2 10
                         (zipmap (keys (first data2))
                                 (keys (first data2)))
                         :YES))


(def forest (read-string (slurp "mRNA_gene_recurrence_forest.edn")))
(show-tree (:tree (last forest)))
