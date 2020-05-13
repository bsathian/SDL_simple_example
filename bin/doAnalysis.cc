#include "doAnalysis.h"

// ./process INPUTFILEPATH OUTPUTFILE [NEVENTS]
int main(int argc, char** argv)
{

//********************************************************************************
//
// 1. Parsing options
//
//********************************************************************************

    // cxxopts is just a tool to parse argc, and argv easily

    // Grand option setting
    cxxopts::Options options("\n  $ doAnalysis",  "\n         **********************\n         *                    *\n         *       Looper       *\n         *                    *\n         **********************\n");

    // Read the options
    options.add_options()
        ("i,input"          , "Comma separated input file list OR if just a directory is provided it will glob all in the directory BUT must end with '/' for the path", cxxopts::value<std::string>())
        ("t,tree"           , "Name of the tree in the root file to open and loop over"                                             , cxxopts::value<std::string>())
        ("o,output"         , "Output file name"                                                                                    , cxxopts::value<std::string>())
        ("n,nevents"        , "N events to loop over"                                                                               , cxxopts::value<int>()->default_value("-1"))
        ("x,event_index"    , "specific event index to process"                                                                     , cxxopts::value<int>()->default_value("-1"))
        ("g,pdg_id"         , "The simhit pdgId match option (default = 13)"                                                        , cxxopts::value<int>()->default_value("13"))
        ("v,verbose"        , "Verbose mode"                                                                                        , cxxopts::value<int>()->default_value("0"))
        ("l,run_ineff_study", "Write debug ntuples 0 = MDs, 1 = SGs, 2 = TLs, 3 = TCs"                                              , cxxopts::value<int>()->default_value("-1"))
       ("d,debug"          , "Run debug job. i.e. overrides output option to 'debug.root' and 'recreate's the file.") 
        ("h,help"           , "Print help")
        ;

    auto result = options.parse(argc, argv);

    // NOTE: When an option was provided (e.g. -i or --input), then the result.count("<option name>") is more than 0
    // Therefore, the option can be parsed easily by asking the condition if (result.count("<option name>");
    // That's how the several options are parsed below

    //_______________________________________________________________________________
    // --help
    if (result.count("help"))
    {
        std::cout << options.help() << std::endl;
        exit(1);
    }

    //_______________________________________________________________________________
    // --input
    if (result.count("input"))
    {
        ana.input_file_list_tstring = result["input"].as<std::string>();
    }
    else
    {
        std::cout << options.help() << std::endl;
        std::cout << "ERROR: Input list is not provided! Check your arguments" << std::endl;
        exit(1);
    }

    //_______________________________________________________________________________
    // --tree
    if (result.count("tree"))
    {
        ana.input_tree_name = result["tree"].as<std::string>();
    }
    else
    {
        std::cout << options.help() << std::endl;
        std::cout << "ERROR: Input tree name is not provided! Check your arguments" << std::endl;
        exit(1);
    }

    //_______________________________________________________________________________
    // --debug
    if (result.count("debug"))
    {
        ana.output_tfile = new TFile("debug.root", "recreate");
    }
    else
    {
        //_______________________________________________________________________________
        // --output
        if (result.count("output"))
        {
            ana.output_tfile = new TFile(result["output"].as<std::string>().c_str(), "create");
            if (not ana.output_tfile->IsOpen())
            {
                std::cout << options.help() << std::endl;
                std::cout << "ERROR: output already exists! provide new output name or delete old file. OUTPUTFILE=" << result["output"].as<std::string>() << std::endl;
                exit(1);
            }
        }
        else
        {
            std::cout << options.help() << std::endl;
            std::cout << "ERROR: Output file name is not provided! Check your arguments" << std::endl;
            exit(1);
        }
    }

    //_______________________________________________________________________________
    // --print_conn
    //_______________________________________________________________________________
    // --nevents
    ana.n_events = result["nevents"].as<int>();
    ana.specific_event_index = result["event_index"].as<int>();

    // -1 upto mini-doublet is all-comb
    // -2 upto segment is all-comb
    // -3 upto tracklet is all-comb NOTE: MEMORY WILL BLOW UP FOR HIGH PU
    // -4 upto trackcandidate is all-comb NOTE: MEMORY WILL BLOW UP FOR HIGH PU
    //  0 nothing
    //  1 upto mini-doublet is all-comb
    //  2 upto mini-doublet is default segment is all-comb
    //  3 upto segment is default tracklet is all-comb
    //  4 upto tracklet is default trackcandidate is all-comb

    //_______________________________________________________________________________
    // --pdg_id
    ana.pdg_id = result["pdg_id"].as<int>();

    //_______________________________________________________________________________
    // --nsplit_jobs

    //_______________________________________________________________________________
    // --verbose
    ana.verbose = result["verbose"].as<int>();

    //
    // Printing out the option settings overview
    //
    std::cout <<  "=========================================================" << std::endl;
    std::cout <<  " Setting of the analysis job based on provided arguments " << std::endl;
    std::cout <<  "---------------------------------------------------------" << std::endl;
    std::cout <<  " ana.input_file_list_tstring: " << ana.input_file_list_tstring <<  std::endl;
    std::cout <<  " ana.output_tfile: " << ana.output_tfile->GetName() <<  std::endl;
    std::cout <<  " ana.n_events: " << ana.n_events <<  std::endl;
    std::cout <<  " ana.run_eff_study: " << ana.run_eff_study <<  std::endl;
    std::cout <<  " ana.run_ineff_study: " << ana.run_ineff_study <<  std::endl;
    std::cout <<  " ana.run_mtv_study: " << ana.run_mtv_study <<  std::endl;
    std::cout <<  " ana.print_centroid: " << ana.print_centroid <<  std::endl;
    std::cout <<  " ana.print_conn: " << ana.print_conn <<  std::endl;
    std::cout <<  " ana.print_boundary: " << ana.print_boundary <<  std::endl;
    std::cout <<  " ana.nsplit_jobs: " << ana.nsplit_jobs <<  std::endl;
    std::cout <<  " ana.job_index: " << ana.job_index <<  std::endl;
    std::cout <<  "=========================================================" << std::endl;

    // Consistency check
    if ((ana.run_ineff_study and ana.run_eff_study) or
        (ana.run_ineff_study and ana.run_mtv_study) or
        (ana.run_eff_study and ana.run_mtv_study)
       )
    {
        RooUtil::error("More than one of -e, -l, or -m option is set! Please only set one of them");
    }


//********************************************************************************
//
// 2. Opening input baby files  (File I/O)
//
//********************************************************************************

    // Create the TChain that holds the TTree's of the baby ntuples
    ana.events_tchain = RooUtil::FileUtil::createTChain(ana.input_tree_name, ana.input_file_list_tstring);

    // Create a Looper object to loop over input files
    // the "www" object is defined in the wwwtree.h/cc
    // This is an instance which helps read variables in the WWW baby TTree
    // It is a giant wrapper that facilitates reading TBranch values.
    // e.g. if there is a TBranch named "lep_pt" which is a std::vector<float> then, one can access the branch via
    //
    //    std::vector<float> lep_pt = www.lep_pt();
    //
    // and no need for "SetBranchAddress" and declaring variable shenanigans necessary
    // This is a standard thing SNT does pretty much every looper we use
    ana.looper.init(ana.events_tchain, &trk, ana.n_events);

//********************************************************************************
//
// Interlude... notes on RooUtil framework
//
//********************************************************************************

    // Set the cutflow object output file
    ana.cutflow.setTFile(ana.output_tfile);
    // ana.cutflow.addCut("CutWeight", UNITY, UNITY);

    // Create ttree to output to the ana.output_tfile
    ana.output_ttree = new TTree("tree", "tree");

    // Create TTreeX instance that will take care of the interface part of TTree
    ana.tx = new RooUtil::TTreeX(ana.output_ttree);

    // Print cut structure
    ana.cutflow.printCuts();

    // pt_boundaries


    // SDL::endcapGeometry.load("scripts/endcap_orientation_data.txt");
    SDL::endcapGeometry.load("/home/users/phchang/public_html/analysis/sdl/TrackLooper_/scripts/endcap_orientation_data_v2.txt"); // centroid values added to the map
    SDL::tiltedGeometry.load("/home/users/phchang/public_html/analysis/sdl/TrackLooper_/scripts/tilted_orientation_data.txt");
    SDL::moduleConnectionMap.load("/home/users/phchang/public_html/analysis/sdl/TrackLooper_/scripts/module_connection_map_data_10_e0_200_100_pt0p8_2p0_400_pt0p8_2p0_nolossers_dxy35cm_endcaplayer2.txt");

    // // Following maps to compute centroid of each modules
    std::map<unsigned int, std::vector<float>> module_xs;
    std::map<unsigned int, std::vector<float>> module_ys;
    std::map<unsigned int, std::vector<float>> module_zs;

    // connection information
    std::ofstream module_connection_log_output;
    if (ana.print_conn)
        module_connection_log_output.open("conn.txt");

    // module boundary information to be written out in case module boundary info is asked to be printed
    std::ofstream module_boundary_output_info;
    if (ana.print_boundary)
        module_boundary_output_info.open("module_boundary.txt");

    // Write the simhits in a given module to the output TTree
    if (ana.print_boundary)
    {
        ana.tx->createBranch<int>("detId");
        ana.tx->createBranch<int>("subdet");
        ana.tx->createBranch<int>("side");
        ana.tx->createBranch<int>("layer");
        ana.tx->createBranch<int>("rod");
        ana.tx->createBranch<int>("module");
        ana.tx->createBranch<int>("ring");
        ana.tx->createBranch<int>("isPS");
        ana.tx->createBranch<int>("isStrip");
        ana.tx->createBranch<vector<float>>("x");
        ana.tx->createBranch<vector<float>>("y");
        ana.tx->createBranch<vector<float>>("z");
    }

    // Looping input file
    while (ana.looper.nextEvent())
    {

        if (ana.specific_event_index >= 0)
        {
            if ((int)ana.looper.getCurrentEventIndex() != ana.specific_event_index)
                continue;
        }

        // If splitting jobs are requested then determine whether to process the event or not based on remainder
        if (result.count("job_index") and result.count("nsplit_jobs"))
        {
            if (ana.looper.getNEventsProcessed() % ana.nsplit_jobs != (unsigned int) ana.job_index)
                continue;
        }

        if (ana.verbose) std::cout <<  " ana.looper.getCurrentEventIndex(): " << ana.looper.getCurrentEventIndex() <<  std::endl;

        // ***************************
        // Testing SDL Implementations
        // ***************************
        //

        // *************************************************
        // Reconstructed hits and formation of mini-doublets
        // *************************************************

        // Main instance that will hold modules, hits, minidoublets, etc. (i.e. main data structure)
        SDL::Event event;

        // Each SDL::Event object in simtrkevents will hold single sim-track related hits
        // It will be a vector of tuple of <sim_track_index, SDL::Event*>.
        std::vector<std::tuple<unsigned int, SDL::Event*>> simtrkevents;

        TStopwatch my_timer;

        // run_eff_study == 0 then run all the reconstruction

        // Adding hits to modules
        for (unsigned int ihit = 0; ihit < trk.ph2_x().size(); ++ihit)
        {

            for (auto& isimhit : trk.ph2_simHitIdx()[ihit])
            {
                int isimtrk = trk.simhit_simTrkIdx()[isimhit];
                if (trk.sim_pt()[isimtrk] < 2.5)
                    continue;
            }


            event.addHitToModule(
                        // a hit
                        SDL::Hit(trk.ph2_x()[ihit], trk.ph2_y()[ihit], trk.ph2_z()[ihit], ihit),
                        // add to module with "detId"
                        trk.ph2_detId()[ihit]
                        );

        }

        float elapsed = 0;

            // ----------------
        if (ana.verbose != 0) std::cout << "Summary of hits" << std::endl;
        if (ana.verbose != 0) std::cout << "# of Hits: " << event.getNumberOfHits() << std::endl;


        // ----------------
        if (ana.verbose != 0) std::cout << "Reco Mini-Doublet start" << std::endl;
        my_timer.Start();
        event.createMiniDoublets();
        float md_elapsed = my_timer.RealTime();
        if (ana.verbose != 0) std::cout << "Reco Mini-doublet processing time: " << md_elapsed << " secs" << std::endl;
        if (ana.verbose != 0) std::cout << "# of Mini-doublets produced: " << event.getNumberOfMiniDoublets() << std::endl;
        if (ana.verbose != 0) std::cout << "# of Mini-doublets considered: " << event.getNumberOfMiniDoubletCandidates() << std::endl;
            // ----------------

            // ----------------
    }
        
    // ********************************************************************************************
    // Perform various studies with reco events and sim-track-matched-reco-hits-based mini-doublets
    // ********************************************************************************************


    // ************************************************
    // Now fill all the histograms booked by each study
    // ************************************************

    // Fill all the histograms


    // Writing output file
    ana.cutflow.saveOutput();

    // Writing ttree output to file
    ana.output_ttree->Write();

    // The below can be sometimes crucial
    delete ana.output_tfile;
}







