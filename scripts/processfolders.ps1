# Before running be sure to run this in PowerShell:
# Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
# When ran, script prepares TCGA downloads from current directory for processing.

# Set the root directory to the current directory
$rootDirectory = Get-Location

# Get all the subfolders in the current directory
$subfolders = Get-ChildItem -Path $rootDirectory -Directory

# Initialize progress counter
$totalCount = $subfolders.Count
$currentCount = 0

# Loop through each subfolder
foreach ($folder in $subfolders) {
    $currentCount++
    $folderName = $folder.Name

    Write-Progress -Activity "Processing folders" -Status "Processing $folderName ($currentCount of $totalCount)" -PercentComplete (($currentCount / $totalCount) * 100)

    # Look for .tsv files directly in the subfolder
    $tsvFiles = Get-ChildItem -Path $folder.FullName -Filter "*.tsv" -File
    foreach ($tsvFile in $tsvFiles) {

        # Move the .tsv file to the root directory
        $destinationPath = Join-Path -Path $rootDirectory -ChildPath $tsvFile.Name
        Move-Item -Path $tsvFile.FullName -Destination $destinationPath -Force
        Write-Output "Moved $($tsvFile.Name) from $folderName to current directory"

        # Modify the TSV file to remove unnecessary rows
        $content = Get-Content -Path $destinationPath
        $rowIndex = 0
        $filteredContent = foreach ($line in $content) {

            # Skip rows 1, 3, 4, 5, and 6 (1-based indexing)
            if ($rowIndex -notin @(0, 2, 3, 4, 5)) {
                $line
            }
            $rowIndex++
        }

        # Write the updated content back to the file
        Set-Content -Path $destinationPath -Value $filteredContent
        Write-Output "Modified $($tsvFile.Name) to remove specified rows"
    }

    # Clean up the current subfolder
    Remove-Item -Path $folder.FullName -Recurse -Force -ErrorAction SilentlyContinue
    Write-Output "Cleaned up folder $folderName"
}

Write-Output "Completed processing $totalCount folders."
Write-Output "All .tsv files extracted, modified, and folders cleaned successfully."
