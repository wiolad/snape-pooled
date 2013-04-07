(*
Copyright (C) 2011-2012 by Emanuele Raineri

date:$Date$
build:$Rev$
*)
open Arg
exception Continue of string
exception Next
let nchr = ref 0
and theta = ref 0.001
and bigd= ref 0.1
and fold = ref ""
and priortype = ref ""
and trust = ref 1.0
and spectrum = ref false
and fast=ref false
and noextremes= ref false
;;
(***************************
		maths 
****************************)
(*binomial table (n choose k) up to n = 20 *)
let binomial_table= 
[|
[|1.|]; 
[|1.; 1.|]; 
[|1.; 2.; 1.|]; 
[|1.; 3.; 3.; 1.|]; 
[|1.; 4.; 6.; 4.; 1.|]; 
[|1.; 5.; 10.; 10.; 5.; 1.|]; 
[|1.; 6.; 15.; 20.; 15.; 6.; 1.|]; 
[|1.; 7.; 21.; 35.; 35.; 21.; 7.; 1.|]; 
[|1.; 8.; 28.; 56.; 70.; 56.; 28.; 8.; 1.|]; 
[|1.; 9.; 36.; 84.; 126.; 126.; 84.; 36.; 9.; 1.|]; 
[|1.; 10.; 45.; 120.; 210.; 252.; 210.; 120.; 45.; 10.; 1.|];
[|1.; 11.; 55.; 165.; 330.; 462.; 462.; 330.; 165.; 55.; 11.; 1.|]; 
[|1.; 12.; 66.; 220.; 495.; 792.; 924.; 792.; 495.; 220.; 66.; 12.; 1.|]; 
[|1.; 13.; 78.; 286.; 715.; 1287.; 1716.; 1716.; 1287.; 715.; 286.; 78.; 13.; 1.|]; 
[|1.; 14.; 91.; 364.; 1001.; 2002.; 3003.; 3432.; 3003.; 2002.; 1001.; 364.; 91.; 14.; 1.|]; 
[|1.; 15.; 105.; 455.; 1365.; 3003.; 5005.; 6435.; 6435.; 5005.; 3003.; 1365.; 455.; 105.; 15.; 1.|]; 
[|1.; 16.; 120.; 560.; 1820.; 4368.; 8008.; 11440.; 12870.; 11440.; 8008.; 4368.; 1820.; 560.; 120.; 16.; 1.|]; 
[|1.; 17.; 136.; 680.; 2380.; 6188.; 12376.; 19448.; 24310.; 24310.; 19448.; 12376.; 6188.; 2380.; 680.; 136.; 17.; 1.|]; 
[|1.; 18.; 153.; 816.; 3060.; 8568.; 18564.; 31824.; 43758.; 48620.; 43758.; 31824.; 18564.; 8568.; 3060.; 816.; 153.; 18.; 1.|]; 
[|1.; 19.; 171.; 969.; 3876.; 11628.; 27132.; 50388.; 75582.; 92378.; 92378.; 75582.; 50388.; 27132.; 11628.; 3876.; 969.; 171.; 19.; 1.|]; 
[|1.; 20.; 190.; 1140.; 4845.; 15504.; 38760.; 77520.; 125970.; 167960.; 184756.; 167960.; 125970.; 77520.; 38760.; 15504.; 4845.; 1140.; 190.; 20.; 1.|]
|]
;;

(* table for pileup parsing *)
let symbols = Hashtbl.create 7;;
Hashtbl.add symbols "A" ( 1 , 7 );;
Hashtbl.add symbols "C" ( 2 , 8 );;
Hashtbl.add symbols "G" ( 3 , 9 );;
Hashtbl.add symbols "T" ( 4 , 10 );;
Hashtbl.add symbols "N" ( 5 , 11 );;
Hashtbl.add symbols "." ( 0 , 6 );;
Hashtbl.add symbols "," ( 0 , 6 );;

(*
rows indexed by k (which can go from 0 to nchr), columns by f 
stores: (n choose k) p^k (1-p)^(n-k)
*)

let prob_k_f = Array.init 1000 (function i -> Array.create 101 0.0)
;; 

let epsilon = sqrt epsilon_float
;;

let isnan (x : float) = x <> x
;;

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
        let res = f n k in
        (Hashtbl.add t2 (n,k) res;
        res) 
;;

let gammaln ( z : float ) =
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
-. tmp +. log(2.5066282746310005*. !ser /. x) 
;;

let factln  =
    memoize (fun (n:int) -> gammaln (float_of_int n +. 1.0)) 
;;

let logbico  = memoize2 
    (fun ( n:int ) ( k:int ) -> factln n -. factln k -. factln ( n - k) )
;;

let logpow ( e:float ) ( base:float ) = 
	if ( e = 0.0 && base = 0.0 ) then 0.0 
	else ( e *. ( log base ) )  
;;

(* binomial probability  (n choose k) p^k q^(n-k) *)

let pbico (n:int) (k:int) (p:float) (q:float) =
if (n<=20) then 
    binomial_table.(n).(k) *. p ** (float_of_int k) *. q ** (float_of_int (n-k))
    else 
    exp (
    logbico n k +. 
    logpow ( float_of_int k ) p +. 
    logpow ( float_of_int n -. float_of_int k ) q
    ) 
;;

let fill_prob_k_f_matrix  n =
(* n number of chromosomes, user-defined *)
    for k = 0 to n do
        for f = 0 to 100 do
            let p = (float_of_int f /. 100.0) in
            prob_k_f.(k).(f) <- ( pbico n k p (1.0 -. p) )  
        done
 done;
;;

let float_of_quality  (c:int) =  
  (10.0**(-.(float_of_int ( c - 33))/. 10.0))
;;

(* 
probability of having 1 alternative allele appearing
in the pileup
    k number of chromosomes bearing an alternative allele
    n total number of chromosomes
    epsilon_list  epsilon_list[0] quality of the ref
                         epsilon_list[1] quality of the alt
*)

let pk ( k:int ) ( n:int ) ( er:float ) ( ea:float )  =
    if (k>n) then  raise (Continue("pk::k>n"));
    let er = if ( er > 0.5 ) then 0.5 else er 
	and ea = if ( ea > 0.5 ) then 0.5 else ea in 
	let pa = ( float_of_int k ) /. ( float_of_int n )
        in 
        let rho    = er *. (1.0 -. 2.0 *. ea)   /. (1.0 -. ea -. er) 
	and alpha  = ea *. (1.0 -. 2.0 *. er)   /. (1.0 -. ea -. er)
        in
        let res = ((1.0 -. rho) *. pa +. alpha *. (1.0 -. pa)) 
        in if (res > 1.0 ) then 
		raise (Continue("internal:invalid pk er:"^(string_of_float er)^" ea:"^(string_of_float ea))) 
		else res
;;

let round ( x:float ) = int_of_float (floor (x +. 0.5))
;;

let first_binomial  (na:int)  (n:int) (g:int)  (er:float) (ea:float)  =
	let probs = Array.init (n+1) (fun i -> 0.0) in 
	for k = 0 to n  do
	let pk =  pk k n er ea  in
	probs.(k)<-( pbico g na pk (1.0 -. pk) );
	done;
	probs
;;

(* 
reader: I know that below I am using a sloppy way of comparing
float numbers. It happens to work in this specific occasion 
*)

let prior_unfolded_informative (theta:float) (beta:float) (bigd:float) (f:float) =
	match f with
	|0.0 -> 1.0 -. theta *. beta -. bigd
	|1.0 -> bigd
	|_   -> theta /. (100.0 *. f)  	 
;;   

let prior_folded_informative theta beta f =
	match f with 
	|0.0 -> ( 1.0 -. theta *. beta ) /. 2.0
	|1.0 -> ( 1.0 -. theta *. beta ) /. 2.0
	|_ -> theta /. (200.0 *. f *. (1.0 -. f))
;;

let prior_unfolded_flat theta bigd f =
	match f with
	|0.0 -> 1.0 -. theta -. bigd
	|1.0 -> bigd
	|_ -> theta /. 99.0
;;

let prior_folded_flat theta f =
	if (f=0.0 || f=1.0) then (1.0 -. theta) /. 2.0
	else  theta /. 99.0 
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
    let lowerbound , upperbound  = 
    if ( !noextremes ) then 1 , ( lp -2 )  
    else 0 , ( lp -1 ) in
    for i = lowerbound to upperbound do
	ef:=!ef +. p.(i) *. f.(i)
    done;
	!ef
;; 

(*****************************************
	parsing
******************************************)

let split_tab=Str.split (Str.regexp "\t") ;;

let fields line =
    let all= split_tab line in
        let get = List.nth all in
    (* chr position ref g pileup qualities *)
    [|get 0;get 1;get 2; get 3; get 4; get 5|]
;;

let parse_pileup pileup qualities =
	let table = Array.init 12 (fun e -> 0) in
	(* table contains the following fields:
	0 number (#) of '.' or ',' characters (reference)
	1 #A
	2 #C
	3 #G
	4 #T	
	5 #N
	6 sum of the quality codes for reference symbols
	7 sum of quality codes for A
        8 sum of quality codes for C
        9 sum of quality codes for G
        10 sum of quality codes for T
        11 sum of quality codes for N
	*)
	let le = String.length pileup in
	if (le = 0) then raise (Continue "empty pileup");
	let i = ref 0 and pos = ref 0 in
        while (!i < le) do
		try (* exception Next *)
		let c = String.uppercase (String.sub pileup !i 1) in
                if (Hashtbl.mem symbols c) then begin
		    let index1 = fst ( Hashtbl.find symbols c )
                    and index2 = snd ( Hashtbl.find symbols c  ) in
                    begin
                        table.(index1)  <- table.(index1) +1 ;
		        table.(index2)  <- table.(index2) + Char.code ( qualities.[!pos] ) ;
		        incr pos ; 
                    end;
                end;
		if (c="$") then raise Next; (* end of read *)
		if (c="^") then begin  incr i; raise Next end; 
                (* ^ followed by quality marks the beginning of a read *)  
		incr i	
		with Next -> begin incr i end;
	done;
	table
;;

let print_spectrum div ps_normalized  =
    let delta = 1.0 /. (float_of_int div) in
  for i = 0 to div-1  do
    Printf.fprintf stdout "%g:%g\t" (delta *. (float_of_int i)) ps_normalized.(i)
  done;
  Printf.fprintf stdout "%g:%g\n" (delta *. (float_of_int div)) ps_normalized.(div) 
;;

let splash () =
  Printf.fprintf stderr "***************************************************************\n%!"; 
  Printf.fprintf stderr "snape-pooled : a method for calling SNPs in pooled samples\n%!";
  Printf.fprintf stderr "$Date$ $Rev$\n%!";
  Printf.fprintf stderr "***************************************************************\n%!"
;;

let parsecmdline () =
 Arg.parse [
    ("-nchr",Set_int (nchr) ,"number of alleles in the pool");
    ("-theta",Set_float (theta) ,"theta");
    ("-D",Set_float (bigd),"D");
    ("-fold",Set_string (fold),"folded or unfolded");
    ("-priortype",Set_string(priortype), "informative or flat");
    ("-trust",Set_float(trust),"[0,1] trust in the reference");
    ("-spectrum",Unit(fun () -> spectrum:= true),"prints MAF probability
    distribution");
    ("-noextremes", Unit(fun () -> noextremes:= true), "excludes f=0 and f=1
    from E(f)")
	] 
    (fun e ->  () ) 
    "Usage e.g. : snape-pooled -nchr 10 -theta 0.001 -D 0.1 -fold folded -priortype informative -spectrum "
;;

let decode_genotype unsorted refc sorted =
        (* calls the genotype looking at the two most frequent nucleotides.
        Badly done, to be replaced by SNAPE model *)
	(* 
	transfer unsorted in an array, as I have to change
	it in place 
	*)
	let genotype = ref "" in
	let copy_unsorted = Array.copy (Array.of_list unsorted) in
	(* first allele *)
	let first = List.nth sorted 0 
	and letters = [| refc;"A";"C";"G";"T" |] in 
	let i = ref 0 in
	try
	while !i <= ( List.length unsorted - 1 ) do
		if ( first > 0 && first = copy_unsorted.(!i) ) then 
		begin
			copy_unsorted.(!i) <- -1; 
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
	while !i <= (List.length unsorted - 1) do
		if ( second >0 && second = copy_unsorted.(!i) ) then 
		begin
			copy_unsorted.(!i) <- -1; 
			raise Exit
		end;
		incr i
	done;
	!genotype
	with Exit -> begin 
            genotype :=!genotype^letters.(!i) 
            end;
	!genotype
;;

let compute_ps_fs div binomial1 prob_k_f prior nchr =
        (* returns normalized ps *)
	let norm=ref 0.0 in
	let fs = Array.init (div +1) (fun e -> 0.0)
		  and ps = Array.init (div +1) (fun e -> 0.0) in
    for j = 0 to div do
			ps.(j)<-0.0;
			fs.(j) <- ( float_of_int j ) /. ( float_of_int div ) ;	
			for i = 0 to nchr do
				ps.(j) <- ps.(j) +. binomial1.(i) *. prob_k_f.(i).(j) ;
			done;
			ps.(j) <- ps.(j) *. (prior fs.(j));
			norm := !norm +. ps.(j) 
		done; 
		for i = 0 to div do
			ps.(i) <- ( ps.(i) /. !norm );
		done;
	[|ps;fs|]
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
      let prior = match fold, priortype with
      |"unfolded","informative" ->  prior_unfolded_informative theta beta bigd  
      |"folded","informative" ->    prior_folded_informative theta beta   
      |"unfolded","flat"        ->  prior_unfolded_flat theta  bigd  
      |"folded","flat"        ->    prior_folded_flat theta   
      |_ -> failwith "unsupported prior" in
    fill_prob_k_f_matrix  nchr;      
    try
       while true do
        let line = input_line inchannel  in
        try 
        let fields = fields line in
        if (fields.(2)="*" ||  (*if the reference base is not known *) 
		(Str.string_match (Str.regexp ".*[+-][0-9].*") fields.(4) 0)  
                (* or there are deletions in the pileup *)
         ) 
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
		if ( ( nr = List.nth sorted_nuc 2 ) && 
                    ( String.length genotype = 2 ) ) 
                then 
                begin
			Printf.fprintf stdout "%s\t%s\t%s\t*\n" chr pos refc;
			raise (Continue "reference could be wrong")
		end;
		let qref = (if (nr>0) then (table.(6) / nr) else 0) 
		and sum_q_non_ref = Array.fold_left (+) 0  ( Array.sub table 7 5 ) 
                in
                let qalt = ( if ( na > 0 ) then ( sum_q_non_ref / na )  else 0 ) in
        if (na > 0 && qalt < 37) then
            begin 
            let skip_msg=
            Printf.sprintf "skipping: %d alt nucleotides with low avg quality:%d" na qalt 
            in
            raise (Continue(skip_msg))
            end;
        let g = (na + nr)  in
        let qalt = if (qalt=0) then qref else qalt in 
        let qref = if (qref=0) then qalt else qref in  
        let er,ea =  (trust *. float_of_quality qref , float_of_quality qalt ) in 
		(* here I fill a vector of length n, containing the values of 
                the first binomial for each possible
		value of k; I do it only once per row as opposed to 100 times *)
		let binomial1 = first_binomial na nchr g er ea 
		in 
	  	let  m = compute_ps_fs div binomial1 prob_k_f prior nchr in
                let ps=m.(0)  and fs=m.(1) in
		(* output fields :
			chr pos ref nr na
			qref qalt genotype 1-p(0) p(1) E(f)
			spectrum (optionally)
		*)	
			let ef = expected_value ps fs 
			in 
		 Printf.fprintf stdout "%s\t%s\t%s\t%d\t%d\t%d\t%d\t%s\t%.4g\t%.4g\t%.4g" 
             	chr pos refc nr na qref qalt genotype
				(1.0 -. ps.(0)) ps.(div) ef;
			if spectrum then begin Printf.fprintf stdout "\t"; print_spectrum div ps end
			else Printf.fprintf stdout "\n" 
          with Continue s -> Printf.fprintf stderr "warning:%s in {%s}\n%!" s line 
			(* raised whenever there's a malformed line *)
        done 
     with End_of_file -> ();;
