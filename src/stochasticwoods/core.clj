(ns stochasticwoods.core
  (:use [clojure.repl]))

(defn- square
  [n]
  (let [nd (double n)]
    (* nd nd)))

(defn- internal-gini
  [m n map_keys]
  (let [nd (double n)]
    (- 1
       (apply +
              (map
               (fn [k] (let [vd (double (get m k 0))]
                        (square (/ vd nd))))
               map_keys)))))

(defn gini-impurity
  [things]
  (let [n (count things)
        freqs (frequencies things)]
    (internal-gini freqs n (keys freqs))))

(defn- combined-gini-score
  [ltemap n1 gtmap n2 map_keys]
  (let [n1d (double n1)
        n2d (double n2)
        ntot (+ n1d n2d)
        s1 (double (internal-gini ltemap n1d map_keys))
        s2 (double (internal-gini gtmap n2d map_keys))]
    (list  (+ (* (/ n1d ntot) s1)
              (* (/ n2d ntot) s2))
           s1 s2)))

(defn- map-add!
  [m k v]
  (assoc! m k (+ v (get m k 0))))

(defn- update-best
  [best ltenew nlte gtnew ngt ltelist gtlist ks split]
  (let [score (combined-gini-score ltenew nlte gtnew ngt ks)]
    (if (< (first score) (first best))
      (concat score (list ltelist gtlist))
      best)))

(set! *warn-on-reflection* true)

(time (decision-tree dd))
(defn best-split-numeric
  [cdata freq_map var]
  (let [ks (keys freq_map)
        get_fn (fn [x] 
                 (.nth ^clojure.lang.PersistentVector
                       (first x) var))
        cdatasort (sort-by (comp get_fn) cdata)]
    (loop [gtgroup (apply list (rest cdatasort))
           ltegroup (list (first cdatasort))
           curr (get_fn (first cdatasort))
           ltefreqs (map-add! (transient {})
                              (second (first gtgroup)) 1)
           nlte 1
           gtfreqs (map-add! (transient freq_map)
                             (second (first gtgroup)) -1)
           ngt (dec (count cdatasort))
           best (list 1 1 1 nil nil nil)]
      (let [gtfirst (first gtgroup)
            dnil (nil? gtfirst)
            d (if dnil nil (get_fn gtfirst))
            update (and (not dnil) (not (== ^double curr ^double d)))
            response (second gtfirst)
            bestnew (if update
                      (update-best best ltefreqs
                                   nlte gtfreqs ngt
                                   ltegroup gtgroup
                                   ks curr) best)]
        (if dnil
          bestnew
          (recur (rest gtgroup)
                 (conj ltegroup gtfirst)
                 d
                 (map-add! ltefreqs response 1)
                 (inc nlte)
                 (map-add! gtfreqs response -1)
                 (dec ngt)
                 bestnew))))))

(defn- cat-freqs
  [cdata var]
  (let [mfn (fn [m [k1 k2]]
              (let [m1 (get m k1 {k2 0})
                    m2 (assoc m1 k2 (inc (get m1 k2 0)))]
                (assoc! m k1 m2)))
        get_fn (fn [x] (let [d (first x)
                            r (second x)]
                        (list
                         (.nth ^clojure.lang.PersistentVector d var)
                         r)))]
    (persistent!
     (reduce mfn (transient {})
             (map (partial get_fn) cdata)))))

(defn- categorical-sort
  [cdata var category]
  (let [get_fn (fn [x] (.nth ^clojure.lang.PersistentVector
                            (first x) var))]
    (loop [data cdata
           ltegroup (list)
           gtgroup (list)]
      (let [x (first data)]
        (if (nil? x)
          (list ltegroup gtgroup)
          (if (= category (get_fn x))
            (recur (rest data) (cons x ltegroup) gtgroup)
            (recur (rest data) ltegroup (cons x gtgroup))))))))

;; (defn- categorical-sort
;;   [cdata var category]
;;   (let [get_fn (fn [x] (.nth ^clojure.lang.PersistentVector
;;                             (first x) var))]
;;     (sort-by #(if (= (get_fn %) category) 0 1) cdata)))

(defn best-split-categorical
  [cdata freq_map var]
  (let [n (count cdata)
        cat_freqs (cat-freqs cdata var)
        categories (keys freq_map)
        f (fn [[k m]]
            (let [dmap (merge-with - freq_map m)
                  n2 (apply + (map second m))]
              (cons k (combined-gini-score m n2 dmap n categories))))
        s (reduce (partial min-key second) (map f cat_freqs))]
    (concat (rest s)
            (categorical-sort cdata var (first s))
            (list (first s)))))

(best-split-categorical [[[:a] :male] [[:b] :female] [[:b] :male]]
                        {:male 2 :female 1}
                        0)

(defn- close-enough?
  [a b]
  (let [ad (double a)
        bd (double b)]
    (< (java.lang.Math/abs (- ad bd)) 1e-10)))

(defn- best-split-all-vars
  [data freqs]
  (apply min-key second
         (map #(cons % (best-split-numeric data freqs %))
              (range (count (first (first data)))))))

(defn infer-categorical
  [cdata var]
  )
(defn decision-tree
  ([data]
   (decision-tree data (gini-impurity (map second data))))
  ([data score]
   (if (zero? score)
     (frequencies (map second data))
     (let [freqs (frequencies (map second data))
           [var gi gil gig
            lte_data gt_data
            split_val] (best-split-all-vars data freqs)]
       (if (or (<= score gi) (close-enough? gi score))
         freqs
         (list (list var split_val)
               (decision-tree lte_data gil)
               (decision-tree gt_data gig)))))))

(def foo (take 10000 (map (fn [[x y]] (if (> (+ (* x x) (* y y)) 8)
                                   [[x y] :female]
                                   [[x y] :male]))
                      (repeatedly #(vector (* 4 (rand)) (* 4 (rand)))))))

(time (decision-tree foo))

(def dd (map list data response))


















(decision-tree [[[1] :a] [[2] :b] [[2] :b]])












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
           (map #(if (categorical? (map (fn [row] (row %)) data))
                   {% :categorical}
                   {% :numeric})
                ks))))

(defn infer-possible-vals
  ([data] (infer-possible-vals data
                               (if (map? (first data))
                                 (keys (first data))
                                 (range (count (first data))))))
  ([data ks]
   (let [f (fn [k] (apply hash-set (map #(% k) data)))]
     (apply merge
            (map (fn [k] {k (f k)}) ks)))))
(infer-possible-vals data [1])
(time (dotimes [n 100] (infer-possible-vals data)))
(defn- break-fns
  ([data ks] (break-fns data ks (infer-datatypes data)))
  ([data ks datatypes]
   (let [datasets (infer-possible-vals data ks)
         get_fn (if (map? (first data)) get nth)]
     (fn [n]
       (let [cmp_fn (if (= :categorical (datatypes n)) = <=)]
         (map (fn [p] (list (fn [x] (cmp_fn (get_fn x n) p)) n p))
              (datasets n)))))))

(time (build-decision-tree data response))

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

(defn- rebuild-split-fn
  [{variable :variable
    split :split
    datatype :datatype}]
  (if (= :categorical datatype)
    (fn [x] (= (get x variable) split))
    (fn [x] (<= (get x variable) split))))

(defn- build-split-fn
  [datatypes labels]
  (fn [[f cat_index split]]
    (let [l (labels cat_index)
          t (datatypes cat_index)
          iscat (= :categorical t)]
      {
       ;; :function (if iscat
       ;;             `(fn [~'x] (= (get ~'x ~cat_index) ~split))
       ;;             `(fn [~'x] (<= (get ~'x ~cat_index) ~split)))
       :label (if iscat
                (str l " == " split)
                (str l " <= " split))
       :variable cat_index
       :split split
       :datatype t})))

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
   (let [close_enough (fn [a b] (< (java.lang.Math/abs (- a b))
                                  1e-5))
         gi0 (double (gini-impurity response))]
     (if (close_enough gi0 0.0)
       {:freqs (frequencies response)
        :gini-impurity gi0}
       (let [m (count data)
             n (count columnkeys)
             bfn (break-fns data columnkeys datatypes)
             breakfns (mapcat bfn columnkeys)
             allbreaks (grab-details data response breakfns)
             splitfn (build-split-fn datatypes labels)
             [mf
              [[lpd lpr lv]
               [rpd rpr rv]]
              gi :as d] (reduce (partial min-key last) allbreaks)]
         (cond
           (close_enough gi0 gi) {:freqs (frequencies response)
                                  :gini-impurity gi0}
           :else (list {:split-fn (splitfn mf)
                        :gini-impurity gi0
                        :freqs (merge-with + lv rv)}
                       (build-decision-tree (vec lpd) (vec lpr)
                                            columnkeys labels datatypes)
                       (build-decision-tree (vec rpd) (vec rpr)
                                            columnkeys labels datatypes))))))))

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
            f (rebuild-split-fn (root :split-fn))]
        (if (f x)
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

;; (defn predict-from-tree
;;   ([tree x] (predict-from-tree tree x eval))
;;   ([tree x eval_fn]
;;    (if (map? tree) (:freqs tree)
;;        (let [[root left right] tree
;;              f (get root :internal_fn)
;;              result (if (nil? f)
;;                       ((eval_fn (get-in root[:split-fn :function])) x)
;;                       (f x))]
;;          (if result
;;            (recur left x eval_fn)
;;            (recur right x eval_fn))))))

;; (defn predict-from-forest
;;   ([forest x]
;;    (predict-from-forest forest x (memoize eval)))
;;   ([forest x eval_fn]
;;    (reduce (partial merge-with +)
;;            (pmap #(predict-from-tree (:tree %) x eval_fn) forest))))

;; (defn classify-from-tree
;;   ([tree x] (classify-from-tree tree x eval))
;;   ([tree x eval_fn]
;;    (first (reduce (partial max-key second)
;;                   (predict-from-tree tree x eval_fn)))))

;; (defn classify-from-forest
;;   [forest x]
;;   (first (reduce (partial max-key second)
;;                  (predict-from-forest forest x))))


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
        ntrees 3000
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
          (pmap #(do (println %2)
                     (hash-map :features (second %1)
                               :indices (first %1)
                               :tree (f %1)))
                (repeatedly samplefn)
                (range)))))

;; (defn- reify-fns-in-tree
;;   ([tree] (reify-fns-in-tree tree eval))
;;   ([tree eval_fn]
;;    (if (map? tree) tree
;;        (let [[root left right] tree
;;              fun (get-in root [:split-fn :function])
;;              f (eval_fn fun)]
;;          (list (assoc root :internal_fn f)
;;                (reify-fns-in-tree left eval_fn)
;;                (reify-fns-in-tree right eval_fn))))))

;; (defn- reify-fns-in-forest
;;   [forest]
;;   (let [meval (memoize eval)
;;         f (fn [tree]
;;             (assoc tree :tree (reify-fns-in-tree (tree :tree) meval)))]
;;     (map f forest)))

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

(defn- shuffle-data
  [data response]
  (let [n (count response)
        s (shuffle (range n))]
    [(mapv (partial get data) s)
     (mapv (partial get response) s)]))

(defn cross-validate
  [data response fold labels category_key]
  (let [[dshuf rshuf] (shuffle-data data response)
        n (count response)
        m (int (/ n fold))
        f (map (fn [x] [(concat (range 0 (* m x))
                               (range (* m (inc x)) n))
                        (range (* m x) (* m (inc x)))])
               (range fold))]
    (map (fn [[train test]]
           (train-and-validate dshuf rshuf train
                               test labels category_key))
         f)))



(defn permute-variable
  [data feature_key]
  (let [v (map #(get % feature_key) data)
        f (fn [row pval] (assoc row feature_key pval))]
    (mapv f data (shuffle v))))

;; (defn- var-imp
;;   [data tree oob_indices]
;;   (let [f (fn [d] (map (fn [i] (classify-from-tree tree (get d i)))
;;                        oob_indices))
;;         c (f data)
;;         f2 (fn [c1 c2] (= c1 c2))]
;;     (fn [feature_key]
;;       (let [m
;;             (apply merge-with
;;                    + (repeatedly 20
;;                                  #(frequencies
;;                                    (map f2 c
;;                                         (f (permute-variable
;;                                             data
;;                                             feature_key))))))]
;;         {feature_key (float
;;                       (/ (get m false 0)
;;                          (+ (get m false 0)
;;                             (get m true 0))))}))))

;; (defn- var-imp-tree
;;   [data]
;;   (let [meval (memoize eval)
;;         all (apply hash-set (range (count data)))]
;;     (fn
;;       [{tree :tree indices :indices features :features}]
;;       (let [t (reify-fns-in-tree tree meval)
;;             ib (apply hash-set indices)
;;             oob (clojure.set/difference all ib)
;;             f (var-imp data t oob)]
;;         (apply merge
;;                (map f features))))))

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
                         (+ (get m false 0)
                            (get m true 0))))}))))

(defn- var-imp-tree
  [data]
  (let [all (apply hash-set (range (count data)))]
    (fn
      [{tree :tree indices :indices features :features}]
      (let [ib (apply hash-set indices)
            oob (clojure.set/difference all ib)
            f (var-imp data tree oob)]
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
    (apply merge-with mean
           (pmap #(do (println %2) (importance_fn %1))
                forest (range)))))

(defn special-read-string
  [x]
  (let [v (read-string x)]
    (if (symbol? v) (keyword v)
        v)))

(def imp_var (time (sort-by second (calc-var-importance data2 forest2))))
(def gini_var (time (sort-by second (calc-var-importance data forest :importance_fn gini-decrease))))

(doseq [a (take 100 (reverse imp_var))]
  ;;(println (gene_list (first a)) (second a))
  (println a)
  )

(doseq [a (sort-by second (time (calc-var-importance data2 forest2 :importance_fn gini-decrease)))]
  (println a))


(let [foo (into [] (sort-by second
                            (time (calc-var-importance data forest
                                                       :importance_fn gini-decrease))))
      n (count foo)]
  (doseq [a (subvec foo (- n 100))]
    ;;(println (gene_list (first a)) (second a))
    (println a)
    ))


(println
 (time (train-and-validate (into [] (take 300 data))
                           (into [] (take 300 response))
                           (into [] (nthrest data 300))
                           (into [] (nthrest response 300))
                           (zipmap (keys (first data))
                                   (keys (first data))))))


(def forest (build-forest data
                          response
                          (zipmap (range) gene_list)
                          ))
(time (count forest))
(println (keys (first data)))



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
(def forest2 (build-forest data2 response2
                           (zipmap (keys (first data2))
                                   (keys (first data2)))))
(count forest2)
(time (build-decision-tree data2 response2))
(def f "/home/kanderson/Code/Bioinformatics/lung_risk/bar2.csv")
(def foo (with-open [rdr (clojure.java.io/reader f)]
           (let [l (line-seq rdr)
                 meta (mapv keyword (take 1000 (clojure.string/split (first l) #", ")))
                 f (fn [row] (mapv special-read-string
                                  (take 1000 (clojure.string/split row #", "))))]
             [meta (doall (mapv f (rest l)))])))
(def foo nil)

(def gene_list (into [] (rest (first foo))))
(def data (map #(into [] (rest %)) (second foo)))
(def response (map first (second foo)))

(def forest (build-forest data
                          response
                          (zipmap (range) gene_list)))
(time (count forest))
(show-tree ((nth forest 4) :tree))


(def foo (cross-validate data2 response2 10 ;;(zipmap (range) gene_list)
                         (zipmap (keys (first data2))(keys (first data2)))
                         :YES))
(def foo2 (cross-validate data response 10 (zipmap (range) gene_list)
                          ;;(zipmap (keys (first data2))(keys (first data2)))
                          :YES))

(calc-auc (nth foo2 3))

(def forest
  (do (nth foo2 4)
      ))
