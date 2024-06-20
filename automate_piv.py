# Developed by John Boamah
import os
import shutil

class PIVDataSorter:
    def __init__(self):
        self.CALIBRATION_FOLDER = "calibration"
        self.RIGHT_FOLDER = "right"
        self.LEFT_FOLDER = "left"
        self.TMP_FILE_EXTENSION = ".tstmp"
        self.LA_FILE_EXTENSION = ".LA.TIF"
        self.LB_FILE_EXTENSION = ".LB.TIF"
        self.RA_FILE_EXTENSION = ".RA.TIF"
        self.RB_FILE_EXTENSION = ".RB.TIF"

    def get_input_folder(self):
        print("Please folder MUST be in current directory")
        main_folder = input("Enter folder name: ")
        return main_folder

    def process_folder(self, main_folder):
        if not main_folder:
            print("Folder name is required. Exiting.")
            return

        print(f'Processing {main_folder}...')

        analysis_folder = os.path.join(main_folder, "Analysis")
        if os.path.exists(analysis_folder):
            shutil.rmtree(analysis_folder)

        run_file = os.path.join(main_folder, f"{main_folder}.run")
        if os.path.exists(run_file):
            os.remove(run_file)

        subfolders = [self.CALIBRATION_FOLDER, self.RIGHT_FOLDER, self.LEFT_FOLDER]
        for subfolder in subfolders:
            os.makedirs(os.path.join(main_folder, subfolder), exist_ok=True)

        raw_data_folder = os.path.join(main_folder, "RawData")
        for filename in os.listdir(raw_data_folder):
            file_path = os.path.join(raw_data_folder, filename)

            if filename.endswith(self.TMP_FILE_EXTENSION):
                os.remove(file_path)
            elif "_Calibration" in filename:
                shutil.move(file_path, os.path.join(main_folder, self.CALIBRATION_FOLDER, filename))
            elif filename.endswith(self.LA_FILE_EXTENSION) or filename.endswith(self.LB_FILE_EXTENSION):
                shutil.move(file_path, os.path.join(main_folder, self.LEFT_FOLDER, filename))
            elif filename.endswith(self.RA_FILE_EXTENSION) or filename.endswith(self.RB_FILE_EXTENSION):
                shutil.move(file_path, os.path.join(main_folder, self.RIGHT_FOLDER, filename))

        print("PIV data sorting completed successfully.")

def main():
    sorter = PIVDataSorter()
    main_folder = sorter.get_input_folder()
    sorter.process_folder(main_folder)

if __name__ == "__main__":
    main()
