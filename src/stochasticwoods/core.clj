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

(defn- break-fns
  ([data ks] (break-fns data ks identity (infer-datatypes data)))
  ([data ks labels] (break-fns data ks labels (infer-datatypes data)))
  ([data ks labels datatypes]
   (let [datasets (infer-possible-vals data ks)]
     (fn [n]
       (map (fn [p]
              (if (= :categorical (datatypes n))
                [(fn [x] (= (get x n) p)) (str "[" (labels n) "] == " p)]
                [(fn [x] (<= (get x n) p)) (str "[" (labels n) "] <= " p)]))
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
  [f n]
  (- 1 (reduce + (map (fn [[k v]] (square (/ v n))) f))))

(defn gini-impurity-after-break
  [[[ld lr lv] [rd rr rv]]]
  (let [nlr (count lr)
        nrr (count rr)
        n (+ nlr nrr)]
    (+ (* (gifromfreqs lv nlr)
          (/ nlr n))
       (* (gifromfreqs rv nrr)
          (/ nrr n)))))

(defn grab-details
  [data response breakfns]
  (map #(let [partition_fn (partition-data-by-breakfn data response)
              p (partition_fn %)]
          (list % p (gini-impurity-after-break p)))
       breakfns))

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
   (let [gi0 (gini-impurity response)
         m (count data)
         n (count columnkeys)
         bfn (break-fns data columnkeys labels datatypes)
         breakfns (mapcat bfn columnkeys)
         allbreaks (grab-details data response breakfns)
         [mf
          [[lpd lpr lv]
           [rpd rpr rv]]
          gi :as d] (reduce (partial min-key last)
                            allbreaks)]
     (cond
       (<= gi0 gi) {:freqs (frequencies response)}
       :else (list {:partition-fn mf :gini-inequality gi
                    :freqs (merge-with + lv rv)}
                   (build-decision-tree (into [] lpd) (into [] lpr)
                                        columnkeys labels datatypes)
                   (build-decision-tree (into [] rpd) (into [] rpr)
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


;; (defn dot-format
;;   ([tree]
;;    (let [al (adjacency-list tree)
;;          edges (map first al)
;;          sides (map #(nth % 2) al)
;;          fs (fn [id node]
;;               (str id
;;                    " [label=\""
;;                    (node :freqs)
;;                    (if (map? node)
;;                      (str node "\"")
;;                      (str (second (node :partition-fn)) "\", shape=\"box\""))
;;                    "];\n"))]
;;      (str
;;       "digraph Tree {\n"
;;       (fs 1 (first tree))
;;       (reduce str
;;               (map (fn [[[e1 e2] node]] (fs e2 node)) al))
;;       (reduce str
;;               (map (fn [[e1 e2] s] (str e1 " -> " e2
;;                                        " [label=\""
;;                                        s
;;                                        "\"];\n"))
;;                    edges sides))
;;       "}"))))

(defn dot-format
  ([tree]
   (let [al (adjacency-list tree)
         edges (map first al)
         sides (map #(nth % 2) al)
         fs (fn [id node]
              (str id
                   " [label=\""
                   (if (contains? node :partition-fn)
                     (str (second (node :partition-fn)) ", "))
                   (node :freqs)
                   (if (contains? node :partition-fn)
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
  [[root & subtree] x]
  (let [nxt {true 0 false 1}
        newtree (nth subtree (nxt ((first (root :partition-fn))  x)))]
    (if (map? newtree) (:freqs newtree) (recur newtree x))))

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
;; (defn train-and-validate
;;   [traindata trainresponse testdata testresponse]
;;   (let [dtree (build-decision-tree traindata trainresponse)
;;         [tpos fpos tneg fneg] [0 0 0 0]
;;         f (fn [x r] (if (= (classify-from-tree dtree x) r)
;;                      :correct
;;                      :incorrect))]
;;     (frequencies (map f testdata testresponse))))

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
                       (+ 2 (int (java.lang.Math/sqrt (count ks)))))
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
          (map #(do (print %2)
                    (hash-map :features (second %1)
                              :indices (first %1)
                              :tree (f %1)))
               (repeatedly samplefn)
               (range)))))

(defn train-and-validate
  [traindata trainresponse testdata testresponse labels]
  (let [forest (build-forest traindata trainresponse labels)
        [tpos fpos tneg fneg] [0 0 0 0]
        f (fn [x r] (if (= (classify-from-forest forest x) r)
                     :correct
                     :incorrect))]
    (frequencies (map f testdata testresponse))))

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

(defn calc-var-importance
  [data forest]
  (let [f (var-imp-tree data)
        mean (fn [a b] (* 0.5 (+ a b)))]
    (apply merge-with mean (map #(do (print %2) (f %1))
                                forest (range)))))

(doseq [a (sort-by second (time (calc-var-importance data forest)))]
  (println a))


(println
 (time (train-and-validate (into [] (take 300 data))
                           (into [] (take 300 response))
                           (into [] (nthrest data 300))
                           (into [] (nthrest response 300))
                           (zipmap (keys (first data))
                                   (keys (first data))))))


(def forest (build-forest data
                          response
                          ;;(zipmap (range) gene_list)
                          (zipmap (keys (first data))
                                  (keys (first data)))
                          ))
(time (count forest))
(println (keys (first data)))

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
                 f (fn [row] (mapv read-string
                                  (clojure.string/split row #", ")))]
             [meta (mapv f (rest l))])))

(def gene_list (into [] (rest (first foo))))
(def data (mapv #(into [] (rest %)) (second foo)))
(def response (mapv first (second foo)))

(def whoa
  (build-decision-tree data response (range (count gene_list))
                       (zipmap (range)
                               (take 300 (nthrest gene_list 300)))))


(def forest (build-forest data
                          response
                          (zipmap (range) gene_list)))
(time (count forest))
(show-tree ((nth forest 1) :tree))

(show-tree (build-decision-tree [[:a 75]
                                 [:a 70]
                                 [:b 55]
                                 [:b 70]]
                                [:true
                                 :false
                                 :false
                                 :true]))
