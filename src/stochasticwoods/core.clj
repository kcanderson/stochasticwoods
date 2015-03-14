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

(defn grab-indices-from-data
  [data indices]
  (let [f (apply juxt (map vector-get-fn indices))]
    (map f data)))

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
     (var-sample-fn (range nfeatures)
                    (int (java.lang.Math/sqrt nfeatures)))))
  ([features mtry]
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
      (sample-without-replacement r (int (/ n 2))))))

(defn sample-fn
  [data]
  (let [n (count data)
        r (range n)]
    (fn []
      (sample-with-replacement
       r n))))

(defn grow-forest
  [data & {:keys [ntrees var_sample_fn decision_tree_fn]
           :or {ntrees 1000
                var_sample_fn (var-sample-fn data)
                decision_tree_fn decision-tree}}]
  (let [score0 (gini-impurity (map second data))
        split_fns (build-split-fns data)
        gi0 (gini-impurity (map second data))
        bootstap_sample_fn (sample-fn data)]
    (take ntrees
          (pmap (fn [inds i]
                 (do (if (zero? (mod i 500))
                       (println i))
                     {:inbag-indices (apply hash-set inds)
                      :tree (decision_tree_fn
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
        f_spec #(/ (get % :true-positpive 0)
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
            (* 0.5
               (+ (- origcorrect
                     (countcorrectfn (permute-var origd variable)))
                  (- origcorrect
                     (countcorrectfn (permute-var origd variable)))))))))))

(defn all-variables-importance
  ([forest data]
   (all-variables-importance
    forest data (range (count (first (first data))))))
  ([forest data variables]
   (let [nt_inv (/ 1.0 (count forest))
         woods (associate-vars-in-forest forest)
         mean (fn [c] (let [n (float (count c))]
                       (if (zero? n)
                         0 (/ (reduce + c) n))))
         vifn (memoize (variable-importance data))
         vivar (fn [v] (list v (* nt_inv
                                 (mean
                                  (remove nil?
                                          (map #((vifn %) v) woods))))))]
     (pmap #(do (if (= 0 (mod %2 5000))
                  (println %2))
                (vivar %1))
           variables (range)))))
;; (defn all-variables-importance
;;   ([forest data]
;;    (all-variables-importance
;;     forest data (range (count (first (first data))))))
;;   ([forest data variables]
;;    (let [woods (associate-vars-in-forest forest)
;;          mean (fn [c] (let [n (float (count c))]
;;                        (if (zero? n)
;;                          0 (/ (reduce + c) n))))
;;          vifn (memoize (variable-importance data))
;;          vivar (fn [v] (list v (mean
;;                                (remove nil?
;;                                        (map #((vifn %) v) woods)))))]
;;      (pmap #(do (if (= 0 (mod %2 5000))
;;                   (println %2))
;;                 (vivar %1))
;;            variables (range)))))
;; (defn all-variables-importance
;;   [forest data]
;;   (let [woods (associate-vars-in-forest forest)
;;         mean (fn [c] (let [n (float (count c))]
;;                       (if (zero? n)
;;                         0 (/ (reduce + c) n))))
;;         vifn (memoize (variable-importance data))
;;         vivar (fn [v] (mean
;;                       (remove nil?
;;                               (map #((vifn %) v) woods))))]
;;     (pmap #(do (if (= 0 (mod % 5000))
;;                  (println %))
;;                (vivar %))
;;           (range (count (first (first data)))))))

;; (defn all-variables-importance
;;   [forest data]
;;   (let [woods (associate-vars-in-forest forest)
;;         mean (fn [c] (let [n (float (count c))]
;;                       (if (zero? n)
;;                         0 (/ (reduce + c) n))))
;;         vifn (memoize (variable-importance data))
;;         vivar (fn [v] (mean
;;                       (remove nil?
;;                               (map #((vifn %) v) woods))))]
;;     (apply concat
;;      (pmap #(map vivar %)
;;            (partition 5000 (range (count (first (first data)))))))))


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


(defn tree-minimum-depth
  ([tree] (tree-minimum-depth tree 0))
  ([tree n]
   (if (leaf? tree) {}
       (merge-with (partial min-key first)
                   {(first (first tree)) (list n)}
                   (tree-minimum-depth (second tree) (inc n))
                   (tree-minimum-depth (third tree) (inc n))))))

(defn- mymean
  [coll]
  (let [n (count coll)]
    (/ (reduce + coll) n)))

(defn minimum-depth
  [forest]
  (apply merge
         (map (fn [[k v]]
                {k (mymean v)})
              (apply merge-with concat
                     (map (comp tree-minimum-depth :tree)
                          forest)))))

(defn root-match?
  [tree m]
  (= (first (first tree)) m))

(defn tree-minimum-depth-interaction
  [tree v1 v2]
  (if (leaf? tree) (list)
      (if (root-match? tree v1)
        (let [lmintree (tree-minimum-depth (second tree))
              rmintree (tree-minimum-depth (third tree))
              f (fn [mt] (if (contains? mt v2)
                          (inc (first (mt v2)))
                          nil))]
          (concat (remove nil? (list (f lmintree) (f rmintree)))
                  (tree-minimum-depth-interaction (second tree) v1 v2)
                  (tree-minimum-depth-interaction (third tree) v1 v2)))
        (concat (tree-minimum-depth-interaction (second tree) v1 v2)
                (tree-minimum-depth-interaction (third tree) v1 v2)))))

(defn minimum-depth-interaction
  [forest v1 v2]
  (map (comp #(tree-minimum-depth-interaction % v1 v2) :tree)
       forest))
