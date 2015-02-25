(ns stochasticwoods.core
  (:use [clojure.repl]))
(set! *warn-on-reflection* true)

(defn sample-with-replacement
  [s n]
  (let [v (if (vector? s) s (into [] s))
        m (count v)]
    (repeatedly n #(nth v (rand-int m)))))

(defn sample-without-replacement
  [things n]
  (take n (shuffle things)))

(defn vector-get-fn
  [var]
  (let [n (int var)]
    (fn [x] (.nth ^clojure.lang.PersistentVector
                 (.first ^clojure.lang.PersistentList x) n))))

(defn third
  [l]
  (second (rest l)))

(defn- square
  [n]
  (let [nd (double n)]
    (* nd nd)))

(defn- internal-gini
  [m n map_keys]
  (let [nd (double n)]
    (- 1
       (apply +
              (map (fn [k] (let [vd (double (get m k 0))]
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

;; (defn- update-best
;;   [best ltenew nlte gtnew ngt ltelist gtlist ks split]
;;   (let [score (combined-gini-score ltenew nlte gtnew ngt ks)
;;         s (double (first score))
;;         b (double (first best))]
;;     (if (< s b(first score)
;;            (first best))
;;       (concat score (list ltelist gtlist split))
;;       best)))

(defn- update-best
  [best ltenew nlte gtnew ngt ltelist gtlist ks split]
  (let [score (combined-gini-score ltenew nlte gtnew ngt ks)
        s (double (first score))
        b (double (first best))]
    (if (< s b)
      (concat score (list ltelist gtlist split))
      best)))

(defn- cmp-fn
  [get_fn]
  (fn [a b]
    (let [da (double (get_fn a))
          db (double (get_fn b))]
      (< da db))))

;; (defn best-split-numeric
;;   [cdata freq_map var]
;;   (let [n (count cdata)
;;         inv (/ 1.0 n)
;;         get_fn (vector-get-fn var)
;;         cdatasort (sort (cmp-fn get_fn) cdata)
;;         r (map second cdatasort)
;;         d (map #(cons %1 %2) (range)
;;                cdatasort)
;;         pd (partition-by (comp get_fn rest) d)
;;         f (fn [i] (let [m (inc i)
;;                        up (take m r)
;;                        down (nthrest r m)]
;;                    (+ (* inv m (gini-impurity up))
;;                       (* inv (- n m) (gini-impurity down)))))
;;         best (apply min-key (comp f first last) pd)
;;         k (inc (first (last best)))
;;         up (take k r)
;;         down (nthrest r k)]
;;     (list (f k)
;;           (gini-impurity up)
;;           (gini-impurity down))
;;     ))

(defn best-split-numeric
  [cdata freq_map var]
  (let [ks (keys freq_map)
        get_fn (vector-get-fn var)
        ;;cdatasort (time (doall (sort (cmp-fn get_fn) cdata)))
        cdatasort (sort-by get_fn cdata)
        ]
    (loop [gtgroup (rest (apply list cdatasort))
           ltegroup (list (first cdatasort))
           curr (double (get_fn (first cdatasort)))
           ltefreqs (map-add! (transient {})
                              (second (first gtgroup)) 1)
           nlte (int 1)
           gtfreqs (map-add! (transient freq_map)
                             (second (first gtgroup)) -1)
           ngt (int (dec (count cdatasort)))
           best (list 1 1 1 nil nil nil)]
      (let [gtfirst (first gtgroup)
            dnil (nil? gtfirst)
            d (if dnil nil (get_fn gtfirst))
            update (and (not dnil) (let [d2 (double d)]
                                 (not (== curr d2))))
            response (second gtfirst)
            bestnew (if update
                      (update-best best ltefreqs
                                   nlte gtfreqs ngt
                                   ltegroup gtgroup
                                   ks curr)
                      best)]
        (cond
          (zero? (first best)) best
          dnil best
          :else (recur (rest gtgroup)
                       (conj ltegroup gtfirst)
                       (double d)
                       (map-add! ltefreqs response 1)
                       (inc nlte)
                       (map-add! gtfreqs response -1)
                       (dec ngt)
                       bestnew))))))

(defn- flatten-persistent!
  [m]
  (let [pm (persistent! m)
        ks (keys pm)]
    (loop [k ks
           m2 (transient pm)]
      (let [k0 (first k)]
        (if (nil? k0) (persistent! m2)
            (recur (rest k)
                   (assoc! m2
                           (first k)
                           (persistent! (get m2 k0)))))))))

(defn- cat-freqs
  [cdata var]
  (let [mfn (fn [m ks]
              (let [k1 (first ks)
                    k2 (second ks)
                    m1 (get m k1 {k2 0})
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

(defn- fast-merge-sub
  [m1 m2]
  (loop [mt1 (transient m1)
         im2 m2]
    (let [i (first im2)]
      (if (nil? i) (persistent! mt1)
          (recur (assoc! mt1 (first i)
                         (- (get mt1 (first i))
                            (second i)))
                 (rest im2))))))


(defn best-split-categorical
  [cdata freq_map var]
  (let [n (count cdata)
        cat_freqs (cat-freqs cdata var)
        categories (keys freq_map)
        f (fn [[k m]]
            (let [dmap (fast-merge-sub freq_map m)
                  n2 (apply + (map second m))]
              (cons k (combined-gini-score m n2 dmap (- n n2)
                                           categories))))
        s (reduce (partial min-key second)
                  (list nil 1 1 1)
                  (map f cat_freqs))]
    (if (== 1 (count cat_freqs)) (list 1 1 1 nil nil nil)
        (concat (rest s)
                (categorical-sort cdata var (first s))
                (list (first s))))))

(defn- close-enough?
  [a b]
  (let [ad (double a)
        bd (double b)]
    (< (java.lang.Math/abs (- ad bd)) 1e-10)))

(defn- min-or-zero
  [data]
  (let [f #(second %)
        z (first (filter #(zero? (f %)) data))]
    (if (nil? z)
      (apply min-key f data)
      z)))

(defn- best-split-all-vars
  ([data freqs split_fns]
   (best-split-all-vars data freqs split_fns
                        (range (count (first (first data))))))
  ([data freqs split_fns vars]
   (min-or-zero (map #(cons % ((nth split_fns %) data freqs %))
                     vars))))

(defn categorical?
  [cdata var]
  (if (every? number? (map (vector-get-fn var) cdata))
    false
    true))

(defn- infer-datatypes
  [cdata]
  (mapv (comp {true :categorical
               false :numeric}
              (partial categorical? cdata))
        (range (count (first (first cdata))))))


;; (defn decision-tree
;;   ([data]
;;    (decision-tree data (range (count (first (first data))))))
;;   ([data vars]
;;    (decision-tree data vars (gini-impurity (map second data))
;;                   (mapv {:categorical best-split-categorical
;;                          :numeric best-split-numeric}
;;                         (infer-datatypes data))))
;;   ([data vars score split_fns]
;;    (if (zero? score)
;;      (frequencies (map second data))
;;      (let [freqs (frequencies (map second data))
;;            [var gi gil gig
;;             lte_data gt_data
;;             split_val] (best-split-all-vars data freqs split_fns vars)]
;;        (if (or (and (<= score gil) (<= score gig) (<= score gi))
;;               (close-enough? gi score))
;;          freqs
;;          (list (list var split_val)
;;                (decision-tree lte_data vars gil split_fns)
;;                (decision-tree gt_data vars gig split_fns)))))))

(defn grab-data-in-order
  [data indices]
  (let [sinds (sort indices)]
    (loop [[i & irest :as inds] sinds
           a 0
           [d & drest :as dtot] data
           ret (list)]
      (if (nil? i) ret
          (if (= a i)
            (recur irest a dtot (cons d ret))
            (recur inds (inc a) drest ret))))))



(defn build-split-fns
  [data]
  (mapv {:categorical best-split-categorical
         :numeric best-split-numeric}
        (infer-datatypes data)))

(defn var-sample-fn
  ([data]
   (let [nfeatures (count (first (first data)))]
     (var-sample-fn data (range nfeatures)
                    (int (java.lang.Math/sqrt nfeatures)))))
  ([data features mtry]
   (fn []
     (sample-without-replacement features mtry))))

(defn decision-tree
  ([data]
   (decision-tree data (var-sample-fn data)))
  ([data vars_sel_fn]
   (decision-tree data vars_sel_fn (gini-impurity (map second data))
                  (build-split-fns data)))
  ([data vars_sel_fn score split_fns]
   (if (zero? score)
     (frequencies (map second data))
     (let [freqs (frequencies (map second data))
           [var gi gil gig
            lte_data gt_data
            split_val] (best-split-all-vars data freqs split_fns
                                            (vars_sel_fn))]
       (if (or (and (<= score gil) (<= score gig) (<= score gi))
              (close-enough? gi score))
         freqs
         (list (list var split_val)
               (decision-tree lte_data vars_sel_fn gil split_fns)
               (decision-tree gt_data vars_sel_fn gig split_fns)))))))

(defn sample-fn
  [data]
  (let [n (count data)
        r (range n)]
    (fn []
      (sample-with-replacement
       (range n) n))))

(defn grow-forest
  [data & {:keys [ntrees var_sample_fn]
           :or {ntrees 1000
                var_sample_fn (var-sample-fn data)}}]
  (let [score0 (gini-impurity (map second data))
        split_fns (build-split-fns data)
        gi0 (gini-impurity (map second data))
        bootstap_sample_fn (sample-fn data)]
    (take ntrees
          (map (fn [inds i]
                 (do (if (zero? (mod i 100))
                       (println i))
                     {:inbag-indices (apply hash-set inds)
                      :tree (decision-tree
                             (grab-data-in-order data inds)
                             var_sample_fn gi0 split_fns)}))
               (repeatedly bootstap_sample_fn) (range)))))

;; (defn grow-forest
;;   [data & {:keys [ntrees features mtry]
;;            :or {ntrees 1000
;;                 features (range (count (first (first data))))}}]
;;   (let [score0 (gini-impurity (map second data))
;;         n (count data)
;;         m (count features)
;;         nfeatures (min m (max (int (+ 5 (java.lang.Math/sqrt m)))
;;                               (int (/ m 10))
;;                               1))
;;         fs (fn [] (list (sample-with-replacement
;;                         (range n) n)
;;                        (sample-without-replacement
;;                         features nfeatures)))
;;         dtypes (infer-datatypes data)
;;         split_fns (mapv {:categorical best-split-categorical
;;                          :numeric best-split-numeric}
;;                         dtypes)
;;         gi0 (gini-impurity (map second data))]
;;     (take ntrees
;;           (map (fn [[inds finds] i]
;;                  (do (println i)
;;                      {:inbag-indices (apply hash-set inds)
;;                       :features finds
;;                       :tree (decision-tree (grab-data-in-order data inds)
;;                                            finds gi0 split_fns)}))
;;                (repeatedly fs) (range)))))

(defn- leaf?
  [item] (or (map? item) (not= 3 (count item))))

(defn split-fn-generator
  ([data] (split-fn-generator data (infer-datatypes data)))
  ([data datatypes]
   (fn [node]
     (let [var (first node)
           split (second node)
           get_fn (vector-get-fn var)]
       (if (= :categorical (.nth ^clojure.lang.PersistentVector
                                 datatypes var))
         (fn [x] (= (get_fn x) split))
         (fn [x] (<= (get_fn x) split)))))))

(defn- freqs-to-prob
  [freqs]
  (let [n (double (apply + (map second freqs)))]
    (apply merge {}
           (map (fn [[k v]]
                  (let [vd (double v)]
                    {k (/ vd n)}))
                freqs))))

(defn resolve-in-tree
  [split_fn_generator]
  (fn [tree x]
    (loop [t tree]
      (if (leaf? t) t
          (if ((split_fn_generator (first t)) x)
            (recur (second t))
            (recur (third t)))))))

(defn classify-in-tree
  [split_fn_generator]
  (fn [tree x]
    (loop [t tree]
      (if (leaf? t) (first (apply max-key second t))
          (if ((split_fn_generator (first t)) x)
            (recur (second t))
            (recur (second (rest t))))))))

(defn oob-votes-tree
  ([data]
   (oob-votes-tree (resolve-in-tree (split-fn-generator data))
                   data))
  ([resolve_fn data]
   (let [indices (apply hash-set (range (count data)))]
     (fn [{tree :tree ib :inbag-indices}]
       (map (fn [i x]
              (if (contains? ib i)
                {}
                (resolve_fn tree x)))
            (range)
            data)))))

(defn oob-votes
  ([forest data]
   (oob-votes (resolve-in-tree (split-fn-generator data))
              forest data))
  ([resolve_fn forest data]
   (let [indices (apply hash-set (range (count data)))
         pf (fn [{tree :tree ib :inbag-indices}]
              (map (fn [i x]
                     (if (contains? ib i)
                       {}
                       (resolve_fn tree x)))
                   (range)
                   data))]
     (reduce #(mapv (partial merge-with +) %1 %2)
             (map pf forest)))))


(defn oob-predictions
  ([forest data]
   (oob-predictions (resolve-in-tree (split-fn-generator data))
                    forest data))
  ([resolve_fn forest data]
   (map freqs-to-prob (oob-votes resolve_fn forest data))))

(defn- fast-merge-add
  [m1 m2]
  (loop [mt1 (transient m1)
         im2 m2]
    (let [i (first im2)]
      (if (nil? i) (persistent! mt1)
          (recur (assoc! mt1 (first i)
                         (+ (get mt1 (first i) 0)
                            (second i)))
                 (rest im2))))))

(defn oob-classifications
  ([forest data]
   (oob-classifications
    (resolve-in-tree (split-fn-generator data)) forest data))
  ([resolve_fn forest data]
   (let [indices (apply hash-set (range (count data)))
         pf (fn [{tree :tree ib :inbag-indices f :features}]
              (map (fn [i x]
                     (if (contains? ib i) {} (resolve_fn tree x)))
                   (range) data))]
     (map ;;(comp first #(apply max-key second %))
      #(if (not (empty? %))
         (first (apply max-key second %))
         nil)
      (reduce #(mapv fast-merge-add %1 %2)
              (map pf forest))))))

(defn classification-performance
  [predictions classifications category threshold]
  (frequencies
   (map #(if (>= (get %1 category 0.5) threshold)
           (if (= %2 category)
             :true-positive
             :false-positive)
           (if (= %2 category)
             :false-negative
             :true-negative)) predictions classifications)))

(defn area-under-curve
  [pts]
  (let [spts (sort-by first pts)
        mean (fn [a b] (* 0.5 (+ a b)))
        area (fn [[[x1 y1] [x2 y2]]]
               (* (mean y1 y2)
                  (- x2 x1)))]
    (apply + (map area (partition 2 1 spts)))))

(defn classification-roc
  [predictions classifications category]
  (let [thresholds (map #(* 0.05 %) (range 20))
        f_spec #(/ (get % :true-positive 0)
                   (+ (get % :true-positive 0)
                      (get % :false-negative 0)))
        fp_rate #(/ (get % :false-positive 0)
                    (+ (get % :false-positive 0)
                       (get % :true-negative 0)))]
    (map (juxt fp_rate f_spec)
         (map (partial classification-performance
                       predictions
                       classifications
                       category)
              thresholds))))


(defn classification-auc
  [predictions classifications category]
  (let [thresholds (map #(* 0.05 %) (range 20))
        f_spec #(/ (get % :true-positive 0)
                   (+ (get % :true-positive 0)
                      (get % :false-negative 0)))
        fp_rate #(/ (get % :false-positive 0)
                    (+ (get % :false-positive 0)
                       (get % :true-negative 0)))]
    (area-under-curve
     (map (juxt fp_rate f_spec)
          (map (partial classification-performance
                        predictions
                        classifications
                        category)
               thresholds)))))

(defn tree-min-var-depth
  ([tree] (tree-min-var-depth tree 1))
  ([tree n]
   (if (leaf? tree) {}
       (let [[[v split] subtree1 subtree2] tree]
         (merge-with #(list (min (first %1) (first %2)))
                     {v (list n)}
                     (tree-min-var-depth subtree1 (inc n))
                     (tree-min-var-depth subtree2 (inc n)))))))

(defn forest-mean-min-depth
  [forest]
  (let [mean (fn [l] (/ (apply + l) (count l)))
        ftree (fn [t] (tree-min-var-depth (t :tree)))]
    (apply merge
           (map #(hash-map (first %) (mean (second %))) 
                (apply merge-with concat
                       (map ftree forest))))))


(defn permute-var
  [data var]
  (let [domain (map (vector-get-fn var) data)
        sdata (shuffle domain)]
    (map #(list (assoc (first %1) var %2) (second %1)) data sdata)))

(defn vars-in-tree
  [tree]
  (if (leaf? tree) (list)
      (cons (first (first tree))
            (concat (vars-in-tree (second tree))
                    (vars-in-tree (third tree))))))

(defn vars-used-in-forest
  [forest]
  (map (comp #(into #{} %) vars-in-tree :tree)
       forest))

(defn associate-vars-in-forest
  [forest]
  (let [vars (vars-used-in-forest forest)]
    (map #(assoc %1 :features %2)
         forest vars)))

(defn variable-importance
  [data]
  (let [orig (into #{} (range (count data)))
        predfn (classify-in-tree (split-fn-generator data))
        classfn (fn [tree d]
                  (map #(= (second %) (predfn tree %)) d))]
    (fn [{tree :tree features :features inbag :inbag-indices}]
      (let [inds (clojure.set/difference orig inbag)
            origd (grab-data-in-order data inds)
            countcorrectfn (fn [d] (count
                                   (filter true? (classfn tree d))))
            origcorrect (countcorrectfn origd)]
        (fn [variable]
          (if (contains? features variable)
            (- origcorrect
               (countcorrectfn (permute-var origd variable)))))))))

(defn all-variables-importance
  [forest data]
  (let [woods (associate-vars-in-forest forest)
        mean (fn [c] (let [n (float (count c))]
                      (if (zero? n)
                        0 (/ (reduce + c) n))))
        vifn (memoize (variable-importance data))
        vivar (fn [v] (mean
                      (remove nil?
                              (map #((vifn %) v) woods))))]
    (map #(do (if (= 0 (mod % 5000))
                (println %))
              (vivar %))
         (range (count (first (first data)))))))


(defn adjacency-list
  ([tree] (adjacency-list tree 1))
  ([tree n]
   (if (leaf? tree)
     [n tree]
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
  ([tree labels datatypes]
   (let [al (adjacency-list tree)
         edges (map first al)
         sides (map #(nth % 2) al)
         fs (fn [id node]
              (str id
                   " [label=\""
                   (if (map? node)
                     (str node)
                     (str (nth labels (first node))
                          (if (= :categorical
                                 (nth datatypes (first node)))
                            (str " == " (second node))
                            (str " <= " (second node)))
                          )) 
                   (if (not (map? node))
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
  [tree labels datatypes]
  (let [dot (dot-format tree labels datatypes)
        pdf (:out (clojure.java.shell/sh
                   "dot" "-Tpdf" :in dot :out-enc :bytes))]
    (clojure.java.shell/sh "display" :in pdf)))


(defn update-dist
  [ntot d]
  (into []
        (concat [0]
                (map (fn [[pp pc] [np nc]]
                       (+ (* pp (/ (- ntot np) ntot))
                          (* pc (/ nc ntot))))
                     (partition 2 1 d)
                     (partition 2 1 (range ntot)))
                (if (< (count d) ntot)
                  [(* (last d) (/ (- ntot (dec (count d))) ntot))]
                  nil))))


(defn special-read-string
  [x]
  (let [v (read-string x)]
    (if (symbol? v) (keyword v)
        v)))



(println (frequencies (map second dd2)))
(println (frequencies (map second dd)))

(def f "/home/kc/Code/Bioinformatics/lung_cancer/data/tumor_rnaseq_normalized.csv")
(def f "/home/kanderson/Code/Bioinformatics/lung_risk/rnaseq_recurrence.csv")
(def f "/home/kanderson/Code/Bioinformatics/lung_risk/lusc_recurrence_rnaseq.csv")
(def f "/home/kanderson/Code/Bioinformatics/lung_risk/luad_metastasis_rnaseq.csv")
(def foo (with-open [rdr (clojure.java.io/reader f)]
           (let [l (line-seq rdr)
                 meta (mapv keyword (clojure.string/split (first l) #", "))
                 f (fn [row] (mapv special-read-string
                                  (clojure.string/split row #", ")))]
             [meta (doall (mapv f (rest l)))])))
(def gene_list (into [] (rest (first foo))))
(def data (map #(into [] (rest %)) (second foo)))
(def response (map first (second foo)))
(def response (map #(cond (= % :M1a) :M1
                          (= % :M1b) :M1
                          :else %)
                   response))
(def foo nil)
(def dd (map list data response))


(def f "/home/kanderson/Code/Bioinformatics/lung_risk/missing_filled_with_year.csv")
(def f "/home/kanderson/Code/Bioinformatics/lung_risk/lusc_clinical_missing_filled.csv")
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
(def features (keys (first data2)))
(def data2 (map #(mapv (fn [k] (get % k))
                       features)
                data2))
(def dd2 (map list data2 response2))

(def f "/home/kanderson/Code/Bioinformatics/lung_risk/luad_clinical_metastasis.csv")
(def f "/home/kanderson/Code/Bioinformatics/lung_risk/lusc_clinical_metastasis.csv")
(def foo (with-open [rdr (clojure.java.io/reader f)]
           (let [l (line-seq rdr)
                 ks (map keyword
                         (clojure.string/split (first l)
                                               #","))
                 f (fn [row] (map special-read-string
                                 (clojure.string/split row #",")))
                 f2 (fn [row] (zipmap ks (f row)))]
             (mapv f2 (rest l)))))
(def datam (mapv #(dissoc % :pathologic_M)
                 foo))
(def responsem (mapv #(cond
                        (= :M1a (:pathologic_M %)) :M1
                        (= :M1b (:pathologic_M %)) :M1
                        :else (:pathologic_M %))
                     foo))
(def featuresm (keys (first datam)))
(def datam (map #(mapv (fn [k] (get % k))
                       featuresm)
                datam))
(def ddm (map list datam responsem))
(def forestm (grow-forest ddm))

(def f "/home/kanderson/Downloads/glop.csv")
(def foo (with-open [rdr (clojure.java.io/reader f)]
           (let [l (line-seq rdr)
                 meta (mapv keyword (clojure.string/split (first l) #", "))
                 f (fn [row] (mapv special-read-string
                                  (clojure.string/split row #", ")))]
             [meta (doall (mapv f (rest l)))])))

(def gene_list (into [] (rest (first foo))))
(def data (map #(into [] (rest %)) (second foo)))
(def response (map first (second foo)))
(def foo nil)
(def dd (map list data response))


(use 'incanter.core)
(use 'incanter.charts)
(use 'incanter.stats)
(view (histogram (map second blah)))
(def x [1 2 5 10 20 30 50 100 200 500 700 900])
(def y (map #(classification-auc
              (oob-predictions (take % woods) dd)
              (map second dd) :YES)
            x))
(def x2 [1 2 5 10 20 50 100 200 500 1000 1500 2000])
(def y2 (map #(classification-auc
               (oob-predictions (take % forest) dd)
               (map second dd) :YES)
             x))
(use 'incanter.pdf)

(view
 (add-lines (xy-plot x y :x-label "Forest Size (num trees)"
                     :y-label "AUC for OOB Recurrence Prediction"
                     :legend true
                     :series-label "mtry = sqrt(p)")
            x2 y2
            :series-label "mtry = p/50"))

(view (histogram
       (let [r (into [] (range 20000))]
         (repeatedly 20000
                     #(count
                       (clojure.set/intersection
                        (into #{} (take 500 (shuffle r)))
                        (into #{} (take 500 (shuffle r)))))))))

(def imp_vars1 (read-string (slurp "rna_important_vars.edn")))
(def imp_vars2 (read-string (slurp "rna_important_vars2.edn")))
(def imp_vars3 (read-string (slurp "rna_important_vars3.edn")))
(def imp_vars4 (read-string (slurp "rna_important_vars4.edn")))
(def imp_vars5 (read-string (slurp "rna_important_vars5.edn")))
(def imp_vars6 (read-string (slurp "rna_important_vars6.edn")))
(def imp_vars7 (read-string (slurp "rna_important_vars7.edn")))
(def imp_vars8 (read-string (slurp "rna_important_vars8.edn")))
(def imp_vars9 (read-string (slurp "rna_important_vars9.edn")))
(def imp_vars10 (read-string (slurp "rna_important_vars10.edn")))

(def woods (grow-forest dd :ntrees 700))
(def imp_vars12 (all-variables-importance woods dd))


(spit "rna_important_vars12.edn" (pr-str imp_vars12))

(def v1 (map (comp #(/ % 5) +)
             imp_vars1 imp_vars3 imp_vars5 imp_vars7 imp_vars9 imp_vars11))
(def v2 (map (comp #(/ % 5) +)
             imp_vars2 imp_vars4 imp_vars6 imp_vars8 imp_vars10 imp_vars12))

(view (scatter-plot imp_vars1 imp_vars2))
(view (histogram (map + v1 v2)
                 :nbins 50))

(defn sig-gene?
  [[x y]]
  (and (> x 0.004) (> y 0.004)))
(defn stringify-gene-keyword
  [k]
  (first (clojure.string/split
          (subs (str k) 1) #"\|")))

(def best_genes
  (map #(stringify-gene-keyword (nth gene_list (last %))) (filter sig-gene? (map list v1 v2 (range)))))
(def best_genes_inds
  (map last (filter sig-gene? (map list v1 v2 (range)))))

(doseq [i best_genes]
  (println i))
(let [chart (scatter-plot v1 v2
                          :x-label "Variable Importance from 6 runs"
                          :y-label "Variable Importance from 6 runs")]
  (view
   (reduce (fn [c [x y i]]
             (if (sig-gene? [x y])
               (add-pointer c x y :text (stringify-gene-keyword
                                         (nth gene_list i)))
               c))
           chart (map list v1 v2 (range)))))

(let [vals (classification-roc
            (oob-predictions redwoods dd)
            (map second dd) :YES)
      vals2 (classification-roc
             (oob-predictions woods dd)
             (map second dd) :YES)]
  (view (add-lines
         (add-lines
          (xy-plot (map first vals2) (map second vals2)
                   :x-label "Specificity"
                   :y-label "Sensitivity"
                   :series-label "Original forest"
                   :legend true)
          (map first vals) (map second vals)
          :series-label "Variable selection then forest")
         [0 1] [0 1] :series-label "Random chance")))





(defn special-sample-fn
  [data var_importance mtry]
  (let [split 0.0025
        n (count var_importance)
        m (count data)
        g1 (map first (filter #(> (second %) split) var_importance))
        n1 (max 1 (int (* mtry (/ (count g1) n))))
        g2 (map first (filter #(<= (second %) split) var_importance))
        n2 (- mtry n1)]
    (fn []
      (list (sample-with-replacement
             (range m) m)
            (concat (sample-without-replacement g1 n1)
                    (sample-without-replacement g2 n2))))))


(def redwoods
  (grow-forest
   dd :samplefn (sample-fn dd best_genes_inds) :ntrees 1000))

(def foo (all-variables-importance redwoods dd))
(map #(nth gene_list (second %))
     (sort-by first > (filter #(> (first %) 0) (map list foo (range)))))
(classification-auc
 (oob-predictions woods dd2)
 (map second dd2) :YES)



(defn scale-map
  [m scale]
  (apply merge
         (map (fn [[k i]] {k (* scale i)})
              m)))
(defn halve
  [n]
  (/ n 2))

(defn count-votes
  [tree]
  (if (leaf? tree) (apply + (map second tree))
      (+ (count-votes (second tree))
         (count-votes (third tree)))))

(defn expand-tree
  [tree var]
  (cond
    (leaf? tree) (list tree)
    (= (first (first tree)) var) (concat
                                  (expand-tree (second tree) var)
                                  (expand-tree (third tree) var))
    :else (for [l (expand-tree (second tree) var)
                r (expand-tree (third tree) var)]
            (list (first tree) l r))))


(defn tree-var-importance
  [data]
  (let [rfn (resolve-in-tree (split-fn-generator data))
        allinds (into #{} (range (count data)))]
    (fn [tree]
      (let [t (tree :tree)
            inbag (tree :inbag-indices)
            oob (clojure.set/difference allinds inbag)
            d (grab-data-in-order data oob)
            preds (fn [t] (map #(rfn t %) d))
            orig (preds t)
            features (tree :features)]
        (reduce merge
                (map
                 (fn [var]
                   (let [e (expand-tree t var)
                         n (count e)
                         votes (reduce #(map (partial merge-with +) %1 %2)
                                       (mapv preds e))]
                     (println (last e))
                     {var (reduce +
                                  (map #(- (to-probability %1 :YES)
                                           (to-probability %2 :YES))
                                       orig votes))}))
                 (take 1 (nthrest features 3))))))))

(defn tree-var-importance
  [data]
  (let [rfn (resolve-in-tree (split-fn-generator data))
        allinds (into #{} (range (count data)))]
    (fn [tree]
      (let [t (tree :tree)
            inbag (tree :inbag-indices)
            oob (clojure.set/difference allinds inbag)
            d (grab-data-in-order data oob)
            preds (fn [t] (map #(rfn t %) d))
            orig (map #(to-probability % :YES) (preds t))
            features (tree :features)]
        (fn [var]
          (let [e (expand-tree t var)
                n (count e)
                votes (reduce #(doall
                                (map (partial merge-with +) %1 %2))
                              (map preds e))]
            ;;(println votes)
            {var (reduce +
                         (map #(- (to-probability %2 :YES)
                                  %1)
                              orig votes))}
            ))))))


(def tree '((:b 1) {:YES 1} ((:a 1) ((:a 1) {:YES 1} {:NO 2}) {:NO 4})))



(defn split-vote
  [split_fn_generator var]
  (letfn [(sv [tree x]
            (loop [t tree]
              (if (leaf? t) t
                  (if (= var (first (first t)))
                    (merge-with (fn [a b] (* 0.5 (+ a b)))
                                (sv (second t) x)
                                (sv (third t) x))
                    (if ((split_fn_generator (first t)) x)
                      (recur (second t))
                      (recur (third t)))))))]
    sv))

(defn to-probability
  [m k]
  (/ (get m k 0)
     (apply + (map second m))))
(defn to-classification
  [m]
  (first (apply max-key second m)))

(defn var-importance
  [data forest]
  (let [orig (map #(to-probability % :YES)
                  (oob-predictions forest data))
        fn_gen (split-fn-generator data)]
    (fn [var]
      (reduce +
              (map #(- %1 (to-probability %2 :YES))
                   orig
                   (oob-predictions 
                    (split-vote fn_gen var)
                    forest data))))))

(defn var-importance
  [data forest]
  (let [orgcls (map second data)
        countfn #(count (filter true? (map = orgcls %)))
        orig (map to-classification
                  (oob-predictions forest data))
        origcorrect (countfn orig)
        fn_gen (split-fn-generator data)]
    (fn [var]
      (- origcorrect
         (countfn (oob-classifications
                   (split-vote fn_gen var)
                   forest data))))))

((var-importance dd2 forest) 1)
(defn all-vars-importance
  [data forest]
  (let [f (var-importance data forest)]
    (map f (range (count (first (first data)))))))

(def forest (grow-forest dd2))
((var-importance dd2 forest) 0)


(def foo (all-vars-importance dd2 forest))
(count foo)

((split-vote (split-fn-generator dd2) 28)
 ((first forest) :tree) (first dd2))


(defn majority-vote
  [votes]
  (first (apply max-key second votes)))

(defn merge-add!
  [tm1 m2]
  (reduce (fn [t [k i]]
            (assoc! t k (+ i (get t k 0))))
          tm1 m2))

(defn maps-reducer
  [maps]
  (map persistent!
       (reduce
        (fn [tm m] (map merge-add! tm m))
        (map transient (first maps))
        (rest maps))))

(defn glop
  [var rfn d forest vars origvotes totorigvotes]
  (reduce (fn [votes update_votes]
            ;;(map #(apply merge-with + %1 %2) votes update_votes)
            votes)
          totorigvotes
          (map (fn [tree tree_vars ov]
                 (if (contains? tree_vars var)
                   (let [newvotes ((oob-votes-tree rfn d) tree)]
                     (map #(if (not= %1 %2)
                             (merge-with - %1 %2))
                          newvotes ov))))
               forest vars origvotes)))
(time (variable-importance dd forest))

(defn variable-importance
  [data forest]
  (let [ntrees (count forest)
        nfeatures (count (first (first data)))
        rfn (resolve-in-tree (split-fn-generator data))
        vars (vars-used-in-forest forest)
        origcls (oob-classifications rfn forest data)
        cls (map second data)
        countcorrectfn (fn [cs]
                         (count (filter true?
                                        (map #(= %1 %2) cs cls))))
        origcount (countcorrectfn origcls)
        origvotes (map #(oob-votes rfn (list %) data) forest)
        totorigvotes (reduce #(map (fn [m1 m2] (merge-with + m1 m2))
                                   %1 %2)
                             origvotes)
        f (fn [d var tree vars_in_tree ov]
            (if (contains? vars_in_tree var)
              (oob-votes rfn (list tree) d)
              ov))
        f2 (fn [d var] (map (partial f d var) forest vars origvotes))
        f3 (fn [var]
             (- origcount
                (countcorrectfn
                 (map majority-vote
                      (maps-reducer
                       (f2 (permute-var data var) var))))))]
    (glop 1000 rfn data forest vars origvotes totorigvotes)))


;; (defn variable-importance
;;   [data forest]
;;   (let [ntrees (count forest)
;;         nfeatures (count (first (first data)))
;;         rfn (resolve-in-tree (split-fn-generator data))
;;         vars (vars-used-in-forest forest)
;;         origcls (oob-classifications rfn forest data)
;;         cls (map second data)
;;         countcorrectfn (fn [cs]
;;                          (count (filter true?
;;                                         (map #(= %1 %2) cs cls))))
;;         origcount (countcorrectfn origcls)
;;         origvotes (map #(oob-votes rfn (list %) data) forest)
;;         f (fn [d var tree vars_in_tree ov]
;;             (if (contains? vars_in_tree var)
;;               (oob-votes rfn (list tree) d)
;;               ov))
;;         f2 (fn [d var] (map (partial f d var) forest vars origvotes))
;;         f3 (fn [var]
;;              (- origcount
;;                 (countcorrectfn
;;                  (map majority-vote
;;                       (maps-reducer
;;                        (f2 (permute-var data var) var))))))]
;;     (map #(do (if (= (mod % 100) 0) (println %)) (f3 %))
;;          (range 20)
;;          ;;(range nfeatures)
;;          )))


(def forest (grow-forest dd :ntrees 2000))
(def imps (variable-importance dd forest))
(first forest)
(def blah (map + (first imps) (second imps)))
(count (second imps))
(def woods (grow-forest dd2 :ntrees 500))
(def imps2 (variable-importance dd2 woods))
(time (count imps2))
(take 20 (sort-by second > (zipmap gene_list imps)))
(println (filter #(> (second %) 0) (zipmap gene_list imps)))
(defn all-variables-importance
  [forest data]
  (let [mean (fn [c] (let [n (count c)]
                      (if (zero? n)
                        0 (/ (reduce + c) n))))
        vifn (variable-importance data)
        vivar (fn [v] (mean
                      (remove nil?
                              (map #(vifn % v) forest))))]
    (map #(do (if (= 0 (mod % 100))
                (println %))
              (vivar %))
         (range (count (first (first data)))))))

(defn vars-in-tree
  [tree]
  (if (leaf? tree)
    #{}
    (let [[[var split] left right] tree
          sub (concat (vars-in-tree left)
                      (vars-in-tree right))]
      (cons var sub))))
(defn unique-vars-in-tree
  [tree]
  (into #{} (vars-in-tree tree)))

(defn shrinking-sample-fn
  [data]
  (let [n (count data)]
    (fn [features]
      (let [m (count features)
            mtry (min m (max 1
                             (java.lang.Math/sqrt m)
                             (/ m 20)
                             ))]
        (list (sample-with-replacement
               (range n) n)
              (sample-without-replacement
               features mtry))))))

(defn increment-map
  [m table]
  (persistent!
   (reduce #(assoc! %1 %2 (inc (get %1 %2 0)))
           (transient m) (into () table))))

(defn remove-entries!
  [m entries]
  (persistent!
   (reduce #(dissoc! %1 %2)
           (transient m) entries)))

(defn grow-forest-shrink-features
  [data & {:keys [ntrees features]
           :or {ntrees 1000
                features (range (count (first (first data))))}}]
  (let [score0 (gini-impurity (map second data))
        split_fns (build-split-fns data)
        gi0 (gini-impurity (map second data))
        samplefn (shrinking-sample-fn data)]
    (loop [n 0
           forest (list)
           fs (into #{} features)
           missing_map {}]
      (if (= n ntrees)
        (list forest fs)
        (let [[inds finds] (samplefn fs)
              sfinds (into #{} finds)
              tree (decision-tree
                    (grab-data-in-order data inds)
                    sfinds gi0 split_fns)
              missing (clojure.set/difference
                       sfinds (unique-vars-in-tree tree))
              m (increment-map missing_map missing)
              r (into #{} (map first (filter #(> (second %) 3) m)))]
          (println n (count fs) (count m)
                   (apply max-key second m))
          (recur (inc n)
                 (cons {:inbag-indices (apply hash-set inds)
                        :features sfinds
                        :tree tree}
                       forest)
                 (clojure.set/difference fs r)
                 m))))))




(def imp_vars (all-variables-importance forestm ddm))
(spit "luad_rnaseq_metastasis.edn" (pr-str imp_vars))
(doseq [a (take 30 (sort-by second > (zipmap gene_list imp_vars)))]
  (println (first a) (second a)))

(println (frequencies (map second ddm)))


(def forest (grow-forest dd))
(def imp_vars (all-variables-importance forest dd))
(count imp_vars)
(classification-auc
 (oob-predictions forest dd)
 (map second dd) :YES)

(defn tree-gini-decrease
  [data]
  (let [all (into #{} (range (count (first (first data)))))]
    (fn [{tree :tree inbag :inbag-indices}]
      (let [oob (clojure.set/difference all inbag)]
        ))))

(defn tree-gini-decrease
  ([data {tree :tree inbag :inbag-indices}]
   (let [all (into #{} (range (count (first (first data)))))
         oob (clojure.set/difference all inbag)]
     (tree-gini-decrease
      data tree (gini-impurity (map second (grab-data-in-order data oob))))))
  ([data {tree :tree inbag :inbag-indices} gi]
   (let [all (into #{} (range (count (first (first data)))))])))


(def forest (grow-forest dd :ntrees 1000))
(def imps (all-variables-importance forest dd))

(def forest12 (grow-forest dd :ntrees 1000
                           :var_sample_fn
                           (var-sample-fn dd finds
                                          (int (java.lang.Math/sqrt(count finds))))))

(def imps12 (all-variables-importance forest12 dd))
(def finds (take 10 (map first
                         (sort-by second > (zipmap (range) imps11)))))
(println (take 20 (sort-by second > (zipmap gene_list imps10))))

(let [z (filter #(> (second %) 0) (last (butlast allimps)))
      z2 (sort-by second > z)]
  (doseq [i (take 20 z2)]
    (println (first i) (/ (second i) 1000))))
(def allimps (map #(zipmap gene_list %) (list imps imps2 imps3 imps4 imps5 imps6 imps7 imps8 imps9 imps10 imps11 imps12)))
(def allforests (list forest forest2 forest3 forest4 forest5 forest6 forest7 forest8 forest9 forest10 forest11 forest12))
(spit "luad_rnaseq_metastasis_impvars_12runs.edn" (pr-str allimps))
(spit "luad_rnaseq_metastasis_forests_12runs.edn" (pr-str allforests))
(def allimps (read-string (slurp "luad_rnaseq_metastasis_impvars_12runs.edn")))
(def allforests (read-string (slurp "luad_rnaseq_metastasis_forests_12runs.edn")))



(defn var-imp-over-runs
  ([all_imps var]
   (map var all_imps)))

(defn var-rank-over-runs
  ([all_imps var]
   (let [imps (var-imp-over-runs all_imps var)]
     (map #(let [foo (zipmap (map first (sort-by second %1)) (range))]
             (if (zero? %2) 0
                 (foo var)))
          all_imps imps))))

(defn grab-top-n-over-runs
  [n all_imps]
  (mapcat #(map first (take n (sort-by second > %))) all_imps))

(defn grab-auc-over-runs
  [data all_forests category]
  (map #(classification-auc
         (oob-predictions % data)
         (map second data) category)
   all_forests))

(def m0classification (grab-auc-over-runs dd allforests :M0))
(def m1classification (grab-auc-over-runs dd allforests :M1))
(def mxclassification (grab-auc-over-runs dd allforests :MX))

(use 'incanter.pdf)
(view
 (add-lines
  (add-lines (xy-plot (range 1 (inc (count allforests)))
                      m0classification
                      :series-label "M0"
                      :x-label "Runs"
                      :y-label "Classification performance (AUC)"
                      :legend true)
             (range 1 (inc (count allforests)))
             m1classification
             :series-label "M1")
  (range 1 (inc (count allforests)))
  mxclassification
  :series-label "MX"))

;; Var imp
(let [mostimps (map first (sort-by second >
                                   (filter #(> (second %) 0)
                                           (last (butlast (butlast (butlast allimps)))))))]
  (let [p (xy-plot (range 1 (inc (count allimps)))
                   (map #(/ % 1000)
                        (var-imp-over-runs allimps (first mostimps)))
                   :series-label (first mostimps)
                   :legend true
                   :x-label "Runs"                   
                   :y-label "VIMP"                   
                   :title "Varible importance over 12 stepwise runs (top 78 genes)")]
    (view (reduce #(add-lines %1
                              (range 1 (inc (count allimps)))
                              (map (fn [i] (/ i 1000))
                                   (var-imp-over-runs allimps %2))
                              :series-label %2)
                  p (rest mostimps)))))

;; Ranks
(let [mostimps (map first
                    (sort-by second >
                             (filter #(> (second %) 0)
                                     (last (butlast (butlast allimps))))))
      ]
  (let [p (xy-plot (range (count allimps))
                   (var-rank-over-runs allimps (first mostimps))
                   :series-label (first mostimps)
                   :legend true
                   :x-label "Runs"                   
                   :y-label "Importance rank"                   
                   :title "Varible rank over 12 stepwise runs (top 38 genes)")]
    (view (reduce #(add-lines %1
                              (range (count allimps))
                              (var-rank-over-runs allimps %2)
                              :series-label %2)
                  p (rest mostimps)))))


(def areas
  (let [genemap (zipmap (range) gene_list)]
    (map #(area-under-curve
           (map list (range)
                (var-imp-over-runs allimps (genemap %))))
     (range (count genemap)))))

(view (histogram areas))
(view (xy-plot (range (count gene_list))
               (map second
                    (sort-by second > (zipmap gene_list areas)))))

(def area_forest
  (let [a (filter #(> (second %) 7.5)
                  (zipmap (range) areas))
        features (map first a)]
    (grow-forest dd :var_sample_fn (var-sample-fn dd features
                                                  (java.lang.Math/sqrt
                                                   (count features))))))

;; (defn collect-stepwise-data
;;   [data]
;;   (loop [vars (range (count (first (first data))))
;;          forests []
;;          imps []
;;          included [(into #{} vars)]]
;;     (if (< (count vars) 10)
;;       (list forests imps included)
;;       (let [n (count vars)
;;             rf (grow-forest
;;                 data
;;                 :var_sample_fn (var-sample-fn
;;                                 data vars (int (java.lang.Math/sqrt n))))
;;             imp (all-variables-importance rf data)
;;             nnext (int (/ n 2))
;;             varsnext (take nnext
;;                            (map first
;;                                 (sort-by
;;                                  second > (zipmap (range) imp))))]
;;         (recur varsnext
;;                (conj forests rf)
;;                (conj imps imp)
;;                (conj included (into #{} varsnext)))))))

(defn collect-stepwise-data
  [data]
  (loop [vars (range (count (first (first data))))
         forests []
         imps []
         included [(into #{} vars)]]
    (let [n (count vars)
          rf (grow-forest
              data
              :var_sample_fn (var-sample-fn
                              data vars (int (java.lang.Math/sqrt n))))
          imp (all-variables-importance rf data)
          varsnext (map first
                        (filter #(> (second %) 0)
                                (zipmap (range) imp)))]
      (println (count vars) (count varsnext))
      (if (= (count vars) (count varsnext))
        (list forests imps included))
      (recur varsnext
             (conj forests rf)
             (conj imps imp)
             (conj included (into #{} varsnext))))))

(def alldata (collect-stepwise-data dd))

(def alldata (collect-stepwise-data dd))
(spit "lusc_rnaseq_recurrence_12runs.edn" (pr-str alldata))

(def alldata (read-string (slurp "lusc_rnaseq_recurrence_12runs.edn")))
(def alldata2 (read-string (slurp "lusc_rnaseq_recurrence_12runs2.edn")))

(def alldata2 (collect-stepwise-data dd))
(spit "lusc_rnaseq_recurrence_12runs2.edn" (pr-str alldata2))
(def allimps (map #(zipmap gene_list %) (second alldata)))
(def allimps2 (map #(zipmap gene_list %) (second alldata2)))

(def allforests2 (first alldata))
(def allimps2 (map #(zipmap gene_list %) (second alldata)))

(def allforests (first alldata))
(def allimps (map #(zipmap gene_list %) (second alldata)))
(spit "lusc_rnaseq_recurrence_impvars_12runs.edn" (pr-str allimps))
(spit "lusc_rnaseq_recurrence_forests_12runs.edn" (pr-str allforests))

(let [mostimps (map first
                    (sort-by second >
                             (filter #(> (second %) 0)
                                     (last (butlast (butlast allimps))))))]
  (let [p (xy-plot (range (count allimps))
                   (var-rank-over-runs allimps (first mostimps))
                   :series-label (first mostimps)
                   :legend true
                   :x-label "Runs"                   
                   :y-label "Importance rank"                   
                   :title "Varible rank over 12 stepwise runs (top 40 genes)")]
    (view (reduce #(add-lines %1
                              (range (count allimps))
                              (var-rank-over-runs allimps %2)
                              :series-label %2)
                  p (rest mostimps)))))

(let [mostimps (map first (sort-by second >
                                   (filter #(> (second %) 0)
                                           (last (butlast (butlast (butlast allimps)))))))]
  (let [p (xy-plot (range 1 (inc (count allimps)))
                   (map #(/ % 1000)
                        (var-imp-over-runs allimps (first mostimps)))
                   :series-label (first mostimps)
                   :legend true
                   :x-label "Runs"                   
                   :y-label "VIMP"                   
                   :title "Varible importance over 12 stepwise runs (top 80 genes)")]
    (view (reduce #(add-lines %1
                              (range 1 (inc (count allimps)))
                              (map (fn [i] (/ i 1000))
                                   (var-imp-over-runs allimps %2))
                              :series-label %2)
                  p (rest mostimps)))))

(def recurrence_classification (grab-auc-over-runs dd allforests :YES))

(view (xy-plot (range 1 (inc (count allforests)))
               recurrence_classification
               :x-label "Runs"
               :y-label "Classification performance (AUC)"
               :title "LUSC recurrence stepwise performance"))


(def impkeys
  (map first
       (filter #(> (second %) 0)
               (last (butlast (butlast allimps))))))

(defn get-best-vimp
  [all_imps var]
  (loop [[imps & rimps] all_imps
         best 0]
    (cond
      (nil? imps) best
      (zero? (get imps var 0)) (recur rimps best)
      :else (recur rimps (get imps var best)))))

(defn mean-vimp
  [all_imps var]
  (mean (map #(get % var) all_imps)))

(defn area-vimp
  [all_imps var]
  (area-under-curve
   (map #(list %1 (get %2 var)) (range) all_imps)))

(defn relative-ranks
  [all_imps vars_used]
  (let [update (fn [imps vars]
                 (let [z (zipmap (range) gene_list)
                       foo (apply merge
                                  (map #(hash-map (z %) (get imps (z %)))
                                       vars))
                       n (count foo)]
                   (reduce (fn [m [[k v] i]]
                             (assoc m k (float (/ i n))))
                           foo
                           (map list
                                (sort-by second > foo)
                                (range)))))
        ranks (map update all_imps vars_used)]
    (fn [var]
      (map #(get % var 0) ranks))))

(defn absolute-ranks
  [all_imps vars_used]
  (let [update (fn [imps vars]
                 (let [z (zipmap (range) gene_list)
                       foo (apply merge
                                  (map #(hash-map (z %) (get imps (z %)))
                                       vars))
                       n (count foo)
                       ;;n2 (- (count imps) n)
                       ]
                   (reduce (fn [m [[k v] i]]
                             (assoc m k i))
                           foo
                           (map list
                                (sort-by second > foo)
                                (range)))))
        ranks (map update all_imps vars_used)]
    (fn [var]
      (map #(get % var 0) ranks))))

(defn relative-rank-area-vimp
  [all_imps vars_used]
  (let [f (relative-ranks all_imps vars_used)]
    (fn [var]
      (area-under-curve
       (map #(list %1 %2) (range) (f var))))))

(defn rank-area-vimp
  [all_imps vars_used rank_fn]
  (let [f (rank_fn all_imps vars_used)]
    (fn [var]
      (area-under-curve
       (map #(list %1 %2) (range) (f var))))))

(def keylist
  (map first
       (take 20
             (sort-by second >
                      (map #(list % (area-vimp allimps %))
                           (keys (first allimps)))))))
(def keylist (keys (first allimps)))

;; Normal RF correlation
(let [k keylist]
  (view
   (scatter-plot
    (map #(/ (get (first allimps) % 0) 1000.0) k)
    (map #(/ (get (first allimps2) % 0) 1000.0) k)
    :x-label "VIMP"
    :y-label "VIMP"
    :title "VIMP correlation for two runs r2=0.0001")))

(let [k keylist]
  (square
   (correlation
    (map #(/ (get (first allimps) % 0) 1000.0) k)
    (map #(/ (get (first allimps2) % 0) 1000.0) k))))


;; "Best" VIMP correlation
(let [k keylist]
  (view
   (scatter-plot
    (map #(get-best-vimp allimps %) k)
    (map #(get-best-vimp allimps2 %) k))))

(let [k keylist]
  (view
   (scatter-plot
    (map #(get-best-vimp allimps %) k)
    (map #(get-best-vimp allimps2 %) k)
    :x-label "VIMP"
    :y-label "VIMP"
    :title "Last VIMP score correlation for two runs r2=0.002")))

;; Mean VIMP score
(let [k keylist]
  (println
   (square
    (incanter.stats/correlation
     (map (partial mean-vimp allimps) k)
     (map (partial mean-vimp allimps2) k)))))

(let [k keylist]
  (view
   (scatter-plot
    (map (partial mean-vimp allimps) k)
    (map (partial mean-vimp allimps2) k)
    :x-label "VIMP"
    :y-label "VIMP"
    :title "Mean VIMP score correlation for two runs r2=0.034")))

;; Area VIMP score
(let [k (keys (first allimps))]
  (println
   (square
    (incanter.stats/correlation
     (map (partial area-vimp allimps) k)
     (map (partial area-vimp allimps2) k)))))

(let [k (keys (first allimps))]
  (println
   (view
    (histogram
     (map (partial area-vimp allimps2) k)
     :nbins 50))))

(let [k keylist]
  (view
   (scatter-plot
    (map (partial area-vimp allimps) k)
    (map (partial area-vimp allimps2) k)
    :x-label "VIMP"
    :y-label "VIMP"
    :title "Area VIMP score correlation for two runs r2=0.057")))

;; Relative rank areas
(let [k keylist
      f3 (rank-area-vimp allimps (third alldata) relative-ranks)
      f4 (rank-area-vimp allimps2 (third alldata2) relative-ranks)]
  (view
   (scatter-plot
    (map f3 k)
    (map f4 k)
    :x-label "VIMP"
    :y-label "VIMP"
    :title "Area VIMP score correlation for two runs r2=0.042")))

;; Absolute ranks
(let [k keylist
      f3 (rank-area-vimp allimps (third alldata) absolute-ranks)
      f4 (rank-area-vimp allimps2 (third alldata2) absolute-ranks)]
  (view
   (scatter-plot
    (map f3 k)
    (map f4 k)
    :x-label "VIMP"
    :y-label "VIMP"
    :title "Area VIMP score correlation for two runs r2=0.37")))

(let [k keylist
      f3 (rank-area-vimp allimps (third alldata) absolute-ranks)      
      f4 (rank-area-vimp allimps2 (third alldata2) absolute-ranks)]
  (square (correlation
           (map f3 k)    
           (map f4 k))))

(let [k keylist
      f3 (rank-area-vimp allimps (third alldata) absolute-ranks)      
      f4 (rank-area-vimp allimps2 (third alldata2) absolute-ranks)]
  (view (histogram
         (map f4 k)
         ;;(map #(* 0.5 (+ %1 %2)) (map f3 k) (map f4 k))
           ;;(map f4 k)
           )))

(let [k (map (zipmap (range) gene_list) (last (butlast (third alldata))))
      f3 (absolute-ranks allimps3 (third alldata))
      ;;f4 (absolute-ranks allimps4 (third alldata2))
      ]
  (let [p (xy-plot (range (count allimps3))
                   (f3 (first k))
                   :series-label (first k)
                   :legend true
                   :x-label "Runs"                   
                   :y-label "Importance rank"                   
                   :title "Varible rank over 12 stepwise runs (top 38 genes)")]
    (view (reduce #(add-lines %1
                              (range (count allimps))
                              (f3 %2)
                              :series-label %2)
                  p (rest k)))))

(doseq [a (sort-by second > (filter #(> (second %) 0) (last (butlast allimps2))))]
  (println a))







