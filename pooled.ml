(*
Copyright (C) 2011 by Emanuele Raineri

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.


*)

open Arg
exception Continue of string
exception Next
(***************************
		maths 
****************************)
let binomial_table= 
[|[|1.|]; [|1.; 1.|]; [|1.; 2.; 1.|]; [|1.; 3.; 3.; 1.|]; [|1.; 4.; 6.; 4.; 
  1.|]; [|1.; 5.; 10.; 10.; 5.; 1.|]; [|1.; 6.; 15.; 20.; 15.; 6.; 
  1.|]; [|1.; 7.; 21.; 35.; 35.; 21.; 7.; 1.|]; [|1.; 8.; 28.; 56.; 70.; 
  56.; 28.; 8.; 1.|]; [|1.; 9.; 36.; 84.; 126.; 126.; 84.; 36.; 9.; 
  1.|]; [|1.; 10.; 45.; 120.; 210.; 252.; 210.; 120.; 45.; 10.; 
  1.|]; [|1.; 11.; 55.; 165.; 330.; 462.; 462.; 330.; 165.; 55.; 11.; 
  1.|]; [|1.; 12.; 66.; 220.; 495.; 792.; 924.; 792.; 495.; 220.; 66.; 
  12.; 1.|]; [|1.; 13.; 78.; 286.; 715.; 1287.; 1716.; 1716.; 1287.; 
  715.; 286.; 78.; 13.; 1.|]; [|1.; 14.; 91.; 364.; 1001.; 2002.; 3003.;
   3432.; 3003.; 2002.; 1001.; 364.; 91.; 14.; 1.|]; [|1.; 15.; 105.; 
  455.; 1365.; 3003.; 5005.; 6435.; 6435.; 5005.; 3003.; 1365.; 455.; 
  105.; 15.; 1.|]; [|1.; 16.; 120.; 560.; 1820.; 4368.; 8008.; 11440.; 
  12870.; 11440.; 8008.; 4368.; 1820.; 560.; 120.; 16.; 1.|]; [|1.; 17.;
   136.; 680.; 2380.; 6188.; 12376.; 19448.; 24310.; 24310.; 19448.; 
  12376.; 6188.; 2380.; 680.; 136.; 17.; 1.|]; [|1.; 18.; 153.; 816.; 
  3060.; 8568.; 18564.; 31824.; 43758.; 48620.; 43758.; 31824.; 
  18564.; 8568.; 3060.; 816.; 153.; 18.; 1.|]; [|1.; 19.; 171.; 969.; 
  3876.; 11628.; 27132.; 50388.; 75582.; 92378.; 92378.; 75582.; 
  50388.; 27132.; 11628.; 3876.; 969.; 171.; 19.; 1.|]; [|1.; 20.; 190.;
   1140.; 4845.; 15504.; 38760.; 77520.; 125970.; 167960.; 184756.; 
  167960.; 125970.; 77520.; 38760.; 15504.; 4845.; 1140.; 190.; 20.; 
  1.|]
|]
;;
let epsilon = sqrt epsilon_float;;
(**
    memoizes any function of one argument
    @param f 'a->'b
    @return a'->'b
*)
let memoize f =
    let t = Hashtbl.create 1000  in 
    fun n ->
		try  Hashtbl.find t n 
    with Not_found ->    
    (* Printf.fprintf stdout  "new2%!"; *) 
        let res = f n  in
        (Hashtbl.add t n res;
        res) 
;; 
let memoize2 f =
    let t2 = Hashtbl.create 1000
    in 
  fun n k ->
    try  Hashtbl.find t2 (n,k) 
    with Not_found ->    
    (* Printf.fprintf stdout  "new2%!"; *) 
        let res = f n k in
        (Hashtbl.add t2 (n,k) res;
        res) 
;; 
let memoize3 f =
    let t3 = Hashtbl.create 1000
    in fun n k p ->
    try  Hashtbl.find t3 (n,k,p)
   with  Not_found -> 
        (* Printf.fprintf stdout  "new3%!"; *)
        let res = f n k p in
        (Hashtbl.add t3 (n,k,p) res;
        res) 
;; 
(**
    logarithm of gamma function
    @param z float
    @return ln |gamma (z)| 
*)
let gammaln z =
    let cof =[|76.18009172947146;-86.50532032941677;
              24.01409824083091;-1.231739572450155;
              0.1208650973866179e-2;-0.5395239384953e-5|]
    in let y = ref z in
    let x = z in
    let tmp=x +.5.5 in
    let tmp = tmp -. (x +. 0.5)*.log(tmp) in
    let ser= ref 1.000000000190015 in
    for j= 0 to 5  do 
        y:=!y +. 1.0;
        ser := !ser +. cof.(j) /. !y 
    done;
-. tmp +. log(2.5066282746310005*. !ser /. x);;
(** ln n! 
    @param n int
    @return ln(n!)
*)
(** ln n!  *)
let factln  =
    memoize (fun n -> gammaln (float_of_int n +. 1.0)) 
;;    
(** 
binomial coefficient n k 
    @param n int
    @param k int
    @return n choose k
*)
let bico n k =
    int_of_float (floor (0.5 +. exp ( factln n -. factln k -. factln ( n - k) )))
;;
let bico = memoize2 bico
;;
let logbico n k = factln n -. factln k -. factln ( n - k)
;;
let logbico  = memoize2 logbico
;;
let logpow e base = if ( e= 0.0 && base = 0.0 ) then 0.0 else e *. log base
;;
(** 
binomial probability distribution
  @param n int
  @param k int
  @param p float
  @param q float
  @return float (n choose k) p^k q^(n-k) 
*)
let pbico n k p q =
	let lnpbico = logbico n k  
  (* ( factln n -. factln k -. factln ( n - k) )*) +. 
  logpow (float_of_int k) p +. logpow (float_of_int n -. float_of_int k) q in
	exp lnpbico;;						
(** 
    from ascii code to epsilon
    @param c int
    @return epsilon float
*)
let float_of_quality  c =  
  (10.0**(-.(float_of_int ( c - 33))/. 10.0))
;;
(** 
probability of having 1 alternative allele appearing
in the pileup
    @param k number of chromosomes bearing an alternative allele
    @param n total number of chromosomes
    @param epsilon_list  epsilon_list[0] quality of the ref
                         epsilon_list[1] quality of the alt
    @return p 
*)
let pk k n er ea  =
    if (k>n) then  raise (Continue("pk::k>n"));
    let er = if (er>0.5) then 0.5 else er 
	and ea = if (ea>0.5) then 0.5 else ea in 
	let pa = (float_of_int k) /. (float_of_int n)
        in 
        let rho1    = er *. (1.0 -. 2.0*.ea)/.(1.0 -. ea -. er) and 
        alpha1  = ea*.(1.0-.2.0*.er)/.(1.0 -. ea -. er)
        in
         let res = ((1.0 -. rho1) *. pa +. alpha1 *. (1.0 -. pa)) 
        in if (res > 1.0) then raise (Continue("internal:invalid pk er:"^(string_of_float er)^" ea:"^(string_of_float ea))) else res ;;      
(** *)
let prob_k n k f =
   if (n<=20) then 
    binomial_table.(n).(k) *. f ** (float_of_int k) *. (1.0 -. f) ** (float_of_int (n-k))   else 
	exp	(logbico n k +. logpow (float_of_int k) f +. logpow (float_of_int (n - k)) (1.0 -. f))  
;;
let prob_k = memoize3 prob_k
;; 
let round x = int_of_float (floor (x +. 0.5))
;;
(** computes p(n_a|f), see docs *)
let p_na_given_f na f n g  er ea =
  let sum = ref 0. in    
    for k = 0 to n  do
    let pk = pk k n er ea in
      sum:=!sum +. (pbico g na pk (1.0 -. pk)) *. (prob_k n k f); 
      (*Printf.fprintf stdout "sum:%g pbico %f prob %f pk %f\n" !sum (pbico g na pk (1.0 -. pk)) (prob_k n k f) pk*)
    done;
  !sum
;;  
(* NB: I know that below I am using a sloppy way of comparing
float numbers. It happens to work in this specific occasion, though*)
let prior_unfolded_informative theta beta bigd f =
	match f with
	|0.0 -> 1.0 -. theta *. beta -. bigd
	|1.0 -> bigd
	|_   -> theta /. (100.0 *. f)  	 
;;   
let prior_folded_informative theta beta f =
	match f with 
	|0.0 -> (1.0 -. theta *. beta ) /. 2.0
	|1.0 -> (1.0 -. theta *. beta ) /. 2.0
	|_ -> theta /. (200.0 *. f *. (1.0 -. f))
;;
let prior_unfolded_flat theta bigd f =
	match f with
	|0.0 -> 1.0 -. theta -. bigd
	|1.0 -> bigd
	|_ -> theta /. 99.0
;;
let prior_folded_flat theta f =
	match f with
	|0.0 -> (1.0 -. theta) /. 2.0 
	|1.0 -> (1.0 -. theta) /. 2.0
	|_ ->  theta /. 99.0 
;;	
(** expected value of f
    @param p float array of probabilities
    @param f float array of frequencies
    @return float
*)
let expected_value p f =
    let lp = Array.length p 
    and lf = Array.length f in
    assert (lp=lf);
	let ef = ref 0.0 in
    for i = 0 to lp -1 do
		ef:=!ef +. p.(i) *. f.(i)
	done;
	!ef
;;                 
(*****************************************
	parsing
******************************************)
let split_tab=Str.split (Str.regexp "\t") ;;
(** 
    @param line
    @return  chr position ref g pileup qualities
*)
let fields line =
    let all= split_tab line in
        let get = List.nth all in
    (* chr position ref g pileup qualities *)
    [|get 0;get 1;get 2; get 3;get 4;get 5|]
let parse_pileup pileup qualities =
	let table = Array.init 8 (fun e -> 0) in
	(* table contains the following fields:
	0 number (#) of '.' or ',' characters (reference)
	1 #A
	2 #C
	3 #G
	4 #T	
	5 #N
	6 sum of the quality codes for reference symbols
	7 sum of quality codes for non reference symbols
	*)
	let le = String.length pileup in
	if (le = 0) then raise (Continue "empty pileup");
	let i = ref 0 and pos = ref 0 in
	while !i<le do
		try
		let c = String.uppercase (String.sub pileup !i 1) in
		if (c="A") then begin 
				table.(1)  <-  table.(1) + 1; 
				table.(7) <- table.(7) + Char.code ( qualities.[!pos] );
				incr pos
				end;
		if (c="C") then begin 
				table.(2)  <-  table.(2) + 1; 
				table.(7) <- table.(7) + Char.code ( qualities.[!pos] );
				incr pos
				end;
		if (c="G") then begin 
				table.(3)  <-  table.(3) + 1; 
				table.(7) <- table.(7) + Char.code ( qualities.[!pos] );
				incr pos
				end;
		if (c="T")  then begin 
				table.(4) <- table.(4) + 1; 
				table.(7) <- table.(7) + Char.code ( qualities.[!pos] );
				incr pos
				end;
		if (c="N")  then begin 
				table.(5) <- table.(5) + 1 ;
				table.(7) <- table.(7) + Char.code ( qualities.[!pos] );
				incr pos
				end;
		if (c="." || c = ",") then begin 
				table.(0) <- table.(0) + 1; 
				table.(6) <- table.(6) + Char.code ( qualities.[!pos] );
				incr pos
				end;
		if (c="$") then raise Next;
		if (c="^") then begin incr i; raise Next end; 
		incr i	
		with Next -> begin incr i end;
	done;
	table
;;
(** *)
let print_spectrum div ps_normalized  =
  for i = 0 to div-1  do
    Printf.fprintf stdout "%g\t" ps_normalized.(i)
  done;
  Printf.fprintf stdout "%g\n" ps_normalized.(div) 
;;
let splash () =
  Printf.fprintf stderr "***************************************************************\n%!"; 
  Printf.fprintf stderr "snape-pooled : a method for calling SNPs in pooled samples\n%!";
  Printf.fprintf stderr "$Date: 2011-08-21 15:50:44 +0200 (Sun, 21 Aug 2011) $ $Rev: 220 $\n%!";
  Printf.fprintf stderr "***************************************************************\n%!"
;;
let nchr = ref 0
and theta = ref 0.001
and bigd= ref 0.1
and fold = ref ""
and priortype = ref ""
and trust = ref 1.0
and spectrum = ref false
;;
let parsecmdline () =
 Arg.parse [
    ("-nchr",Set_int (nchr) ,"number of alleles in the pool");
    ("-theta",Set_float (theta) ,"theta");
    ("-D",Set_float (bigd),"D");
    ("-fold",Set_string (fold),"folded or unfolded");
    ("-priortype",Set_string(priortype), "informative or flat");
	("-trust",Set_float(trust),"[0,1] trust in the reference");
    ("-spectrum",Unit(fun () -> spectrum:= true),"prints MAF probability distribution")
	] 
    (fun e ->  () ) 
    "Usage e.g. : snape-pooled -nchr 10 -theta 0.001 -D 0.1 -fold folded -priortype informative -spectrum "
;;
let decode_genotype unsorted refc sorted =
	(* 
	transfer unsorted in an array, as I have to change
	it in place 
	*)
	let genotype = ref "" in
	let copy_unsorted = Array.copy (Array.of_list unsorted) in
	(* first allele *)
	let first = List.nth sorted 0 
	and letters = [|refc;"A";"C";"G";"T"|] in 
	let i = ref 0 in
	try
	while !i<= (List.length unsorted - 1) do
		if (first >0 && first = copy_unsorted.(!i)) then 
		begin
			copy_unsorted.(!i)<- -1; 
			raise Exit
		end;
		incr i
	done;
	!genotype
	with Exit -> begin genotype :=!genotype^letters.(!i) end;
	i:=0;
	(* second allele *)
	let second =  List.nth sorted 1 in
	if (second = 0) then !genotype else  
	try
	while !i<= (List.length unsorted - 1) do
		if (second >0 && second = copy_unsorted.(!i)) then 
		begin
			copy_unsorted.(!i)<- -1; 
			raise Exit
		end;
		incr i
	done;
	!genotype
	with Exit -> begin genotype :=!genotype^letters.(!i) end;
	!genotype
;;
(** main *)
let _ =
    	splash ();
		parsecmdline ();  
       	let nchr = !nchr
		and theta = !theta
		and bigd = !bigd
		and fold = !fold
		and priortype = !priortype
		and trust = !trust 
		and spectrum = !spectrum  in
      	let bigd = if ( bigd >= theta ) then bigd else theta  in 
      	let inchannel = stdin in  
		  let div=100 in let beta = log (float_of_int div) +. 0.57721 in
		  let fs = Array.init (div +1) (fun e -> 0.0)
		  and ps = Array.init (div +1) (fun e -> 0.0) in
      let prior = match fold, priortype with
      |"unfolded","informative" ->  prior_unfolded_informative theta beta bigd  
      |"folded","informative" ->    prior_folded_informative theta beta   
      |"unfolded","flat"        ->  prior_unfolded_flat theta  bigd  
      |"folded","flat"        ->    prior_folded_flat theta   
      |_ -> failwith "unsupported prior" in
      try
       while true do
        let line = input_line inchannel  in
        try 
        let fields = fields line in
        if (fields.(2)="*" ||  (*if the reference base is not known *) 
		(Str.string_match (Str.regexp ".*[+-][0-9].*") fields.(4) 0)) (* or there are deletions in the pileup *)
		then raise (Continue("skipping: contains indels or unknown ref")); (* skip the line *)
        let chr = fields.(0) and pos = fields.(1)
		and refc = fields.(2) in   
		let table = parse_pileup fields.(4) fields.(5) in
		let na = table.(1) + table.(2) + table.(3) + table.(4) + table.(5) 
		and nr = table.(0)
		in
		let nuc=[table.(0);table.(1);table.(2);table.(3);table.(4)]
		in let sorted_nuc = List.sort (fun x y -> - compare x y) nuc in
		let genotype = decode_genotype nuc refc sorted_nuc
		in 
		if ((nr = List.nth sorted_nuc 2) && (String.length genotype =2 ) ) then begin
			Printf.fprintf stdout "%s\t%s\t%s\t*\n" chr pos refc;
			raise (Continue "reference could be wrong")
		end;
		let qref = (if (nr>0) then (table.(6) / nr) else 0) 
		and qalt = (if (na>0) then (table.(7) / na)  else 0) in 
        let g = (na + nr)  in
        let qalt = if (qalt=0) then qref else qalt in 
        let qref = if (qref=0) then qalt else qref in  
        let er,ea =  (trust *. float_of_quality qref , float_of_quality qalt ) in 
		let norm = ref 0.0 in
		for i = 0 to div do
			fs.(i) <- (float_of_int i) /. (float_of_int div) ;
			ps.(i) <- ((p_na_given_f na fs.(i) nchr g er ea) *. (prior fs.(i)));  
			norm := !norm +. ps.(i)  		
		done;
		for i = 0 to div do
			ps.(i) <- ( ps.(i) /. !norm )
		done;
		(* output fields :
			chr pos ref nr na
			qref qalt genotype 1-p(0) p(1) E(f)
			spectrum (optionally)
		*)	
			let ef = expected_value ps fs 
			in 
		 Printf.fprintf stdout "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%.4g\t%.4g\t%.4g" 
             	chr pos refc nr na qref qalt genotype
				(1.0 -. ps.(0)) ps.(1) ef;
			if spectrum then begin Printf.fprintf stdout "\t"; print_spectrum div ps end
			else Printf.fprintf stdout "\n" 
          with Continue s -> Printf.fprintf stderr "warning:%s in {%s}\n%!" s line 
			(* raised whenever there's a malformed line *)
        done 
     with End_of_file -> ();;
(** *)
