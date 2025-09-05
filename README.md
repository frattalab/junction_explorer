So you want to make some plots. And the splice junctions we hosted at Shiny.IO didn't cut it. 

If we could have freely made all the junctions available, we would have, but the full NYGC dataset is something like 5GB and the RiMOD over a gig as well, and free things have limits. 

But do not fear! You are an intrepid traveller and you're willing to try and deploy this app on your desktop

First, you will go to Zenodo and get the full datasets. 
You will go here 
https://zenodo.org/records/15977268

And download the data here
full_rimod_parquets.zip
nygc_parquets.zip


You are going to download this GitHub repo (this one, where we are, right now)
<img width="930" height="471" alt="Screenshot 2025-09-05 at 12 47 42" src="https://github.com/user-attachments/assets/075cbd0b-b078-4bad-a074-8cd09d58dac8" />


You are going to upzip those folders.

You are going to navigate to 
data > parquets
<img width="1324" height="184" alt="Screenshot 2025-09-03 at 20 08 18" src="https://github.com/user-attachments/assets/8db69cf7-e9cf-4448-885e-3af91a4fe19f" />

You are going to delete what is inside the parquets folder, not the folder itself. 
<img width="328" height="163" alt="Screenshot 2025-09-03 at 20 09 59" src="https://github.com/user-attachments/assets/c381dbc7-9a80-45d1-863b-57c11937b93c" />

You are going to put the thing you have just downloaded (after unzipping) into the the parquets folders. 

You are going to double click on the project
<img width="1042" height="228" alt="Screenshot 2025-09-03 at 20 11 30" src="https://github.com/user-attachments/assets/a528f444-a7f0-42f9-ac15-ab384c8a6b9d" />

RStudio is going to open, and you are not afraid, you can do this. 
You are going to click on the `app.R` to open it
<img width="487" height="400" alt="Screenshot 2025-09-03 at 20 12 46" src="https://github.com/user-attachments/assets/6d68740a-751c-40a3-bea8-d3bcd330abc0" />

You are going to click the button saying run app at the upper right hand side 
<img width="1014" height="475" alt="Screenshot 2025-09-05 at 12 08 29" src="https://github.com/user-attachments/assets/d6347b0a-5e3c-441d-b9c8-1e8cb0a687a3" />



Test on, say, the annotated version of the STMN2 cryptic exon
`chr8:79611214-79636802` which won't show up in the hosted version
